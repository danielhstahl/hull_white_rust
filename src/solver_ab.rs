//! A/B cost comparison: the bracketed rootfinder vs the `nrfind` Newton iteration it replaced.
//!
//! Accuracy of the new solve was checked against quadrature when it landed; cost was not.  This
//! harness runs both solvers against *the same* objective and derivative closures the production
//! Jamshidian code builds (`coupon_bond_generic_t` over `bond_price_t_raw`, and
//! `coupon_bond_price_t_deriv`), over a grid of curve fixtures x coupon schedules x strikes, so
//! the only thing that varies is the solver.
//!
//! What is compared per case:
//!
//! * `f_calls` / `df_calls` — the solver's asks on the model.  A price costs what an evaluation
//!   costs, so this is the machine-independent measure of the solve.
//! * `curve_calls` — how many times the solve reached for the yield/forward curves, counted by
//!   the curves themselves.  Independent of solver internals entirely.
//! * `setup` — curve calls spent building the analytic bracket plus the seed.  The new solver has
//!   a setup step the old one never paid, so it is counted rather than hidden.
//! * `ns` — wall clock per solve, median of repeats, old and new interleaved in the same repeat
//!   loop so clock drift and cache state hit both alike.
//! * root agreement / residual — whether the two solvers answered the same question.
//! * `d_ratio` — the analytic derivative over a central difference of the objective at the root.
//!   If the analytic derivative is not the objective's derivative, the min-progress guard rejects
//!   every Newton step and the solve degenerates into pure bisection.
//!
//! Run it (it is env-gated so the ordinary test suite stays fast):
//!
//! ```text
//! SOLVER_AB=1 cargo test --release solver_ab -- --nocapture   // sections A-D, the two production shapes
//! SOLVER_AB=1 cargo test --release exit_rule_variants -- --nocapture  // section E, candidate exit rules
//! SOLVER_AB=1 cargo test --release old_solver_stress -- --nocapture   // section F, where the old solve breaks
//! AB_TRACE=1  cargo test --release trace_one_solve -- --nocapture     // iteration-by-iteration dump
//! AB_REPS=51  ...                                                  // more repeats for timing
//! ```
//!
//! ## Findings (84 cases: 4 curve shapes x {4,8,48} coupons x 7 strikes, release build)
//!
//! | solver | iters | evals | curve calls | ns/solve | max root err vs 1e-15 ref |
//! |---|---|---|---|---|---|
//! | `nrfind` (what we had) | 4.6 | 832 | 59,688 | ~6,700 | ~9.5e-13 vs new |
//! | `rootfinder::solve` (pre-port) | 53.1 | 9,092 | 575,484 | 92,712 | -- |
//! | `v2` -- exit on the Newton correction | 25.4 | 4,565 | 289,032 | 45,600 | 5.6e-17 |
//! | `v3` -- exit on the price residual | 25.5 | 4,458 | 296,604 | n/a | **9.6e-9 -- unsafe** |
//! | `v4` -- residual-decrease step acceptance, rate-space exits | 5.2 | 1,463 | 101,208 | 15,153 | 2.2e-16 |
//! | `rootfinder::solve` (port of `v4`, what ships now) | 5.4 | 1,442 | 102,132 | ~15,900 | 2.2e-16 |
//!
//! The shipped solve costs ~11x the evaluations it replaced, at the same accuracy.  `v4` lands at
//! ~1.8x `nrfind` (1,463 vs 832 evals, 101,208 vs 59,688 curve calls) -- the extra is the
//! residual probe that buys step acceptance -- while keeping the bracket, the seed and every
//! failure mode the rewrite was done for.  It is not a bad derivative
//! (`d_ratio` = 1.0000 against a central difference) and not a huge bracket
//! (0.13 - 2.6 rate units on this grid).  `AB_TRACE=1` shows where the passes go:
//!
//! * Iteration 1: the Newton step from `mu_r` lands ~2e-4 from the root.
//! * The min-progress guard rejects it -- `|correction| (1.4e-2) < 0.1 x width (3.5e-1)` -- so
//!   the solver bisects instead.
//! * From there it alternates: a near-root step whose correction is tiny (already converged) gets
//!   rejected again, and bisection spends the remaining ~35 passes grinding the *stale far end*
//!   of the bracket down to the 1e-12 rate tolerance.
//!
//! The tolerance is on the rate, which is the right thing to ask for here, so the fix is not to
//! loosen it -- it is to stop rejecting the step that already converged.
//!
//! * `v2` (accept "correction below tolerance" as convergence) recovers ~2x and is exactly
//!   conservative otherwise.
//! * `v3` shows why the exit must stay in rate space: `f` is a price and can be steep, so a
//!   1e-8 price residual hid a 7e-3 error in the critical rate, which would put every leg
//!   strike in the Jamshidian decomposition wrong while the objective reported "zero".
//! * `v4` (accept a Newton step that lands inside the bracket and strictly reduces the residual,
//!   in addition to the width rule; all *exits* still in rate space) matches the old solve's
//!   iteration count -- ~5 -- while keeping the bracket, the seed and every failure mode.
//!
//! `v4`'s cost is bounded where the guard was protecting us: on the deep exponential-tail
//! fixture it takes 119 iterations to the exact root versus the pre-port solver's 65 (bisection's
//! own floor on that fixture is ~48), so the pathological-tail penalty is ~2x, and it is paid
//! only in a regime that does not occur on the tested model grid.  The shipped `max_iterations`
//! still bounds it.  The port removed that penalty rather than accepting it -- see below.
//!
//! ## Ported: what `rootfinder::solve` does now
//!
//! The `v4` rule is in the shipped loop: `inside && (improves || makes_progress)`, plus the
//! correction-based convergence exit from `v2`.  Both exits stay in rate space; the residual is
//! used as a *step acceptance* test and never as an exit, because `v3` showed a small price
//! residual can hide a large error in the rate.
//!
//! `improves` ships as `|f(x_newton)| <= MAX_MODEL_ERROR_RATIO * |f(x)|` with the ratio at
//! `0.25`, not `v4`'s bare "any strict decrease".  That is deliberate and it is the tail that
//! decides it: a locally exponential objective contracts by exactly `exp(-1) = 0.368` per
//! accepted Newton step *no matter how far away the root is*, so a threshold above `1/e` still
//! admits the crawl the width rule was written to hand back to bisection.  Threshold scan
//! (84-case grid; deep-tail fixture `exp(-x) = 1e-6` seeded at -100, release build):
//!
//! | acceptance threshold for `improves` | grid iters avg | grid evals | deep-tail iters |
//! |---|---|---|---|
//! | width rule only (pre-port) | 53.1 | 9,092 | 65 |
//! | any decrease (`v4`) | 5.2 | 1,463 | 119 |
//! | `<= 0.40` | 5.2 | 1,463 | 119 (blows a 100 cap) |
//! | `<= 0.33` | 5.2 | 1,463 | 11 |
//! | **`<= 0.25` (shipped)** | 5.4 | 1,442 | 12 |
//! | `<= 0.10` | 6.6 | 1,720 | 12 |
//! | `<= 0.05` | 7.0 | 1,804 | 14 |
//!
//! Root agreement with a 1e-15 reference is flat across the entire scan (worst `2.3e-16`), so the
//! tighter threshold buys tail safety for ~4% of the solve's evaluations and costs nothing in
//! accuracy.  `0.25` sits a factor 1.47 below `1/e`, clear of the contraction rate rather than
//! sitting on it.
//!
//! After the port (84 grid cases, release): **5.4 iters / 1,442 evals / 102,132 curve calls /
//! ~16us per solve**, against 53.1 / 9,092 / 575,484 / ~93us pre-port and 5.0 / 832 / 59,688
//! for `nrfind` -- 1.73x the old evaluation count instead of 10.9x, with the analytic bracket,
//! the `mu_r` seed and every failure mode intact.  Worst root error against the 1e-15 reference:
//! `2.2e-16` on the grid, `1.8e-14` on the adversarial grid.  Deep tail: 12 iterations, so the
//! tail penalty `v4` carried is gone.
//!
//! Three always-on guards (no env var) hold the port in place:
//!
//! * `ported_step_rule_matches_a_tight_reference_and_stays_cheap` -- root agreement with
//!   `solve_v0@1e-15` within `1e-12`, mean iterations `<= 10` (measured 5.4 vs 53.1 pre-port),
//!   total evaluations within 2x the old `nrfind` solve (measured 1.73x).  The evaluation bound
//!   is exact-count based, so it is machine independent.
//! * `ported_step_rule_holds_on_the_adversarial_grid` -- same root agreement over the 48
//!   K = 1e-6..1e6 cases, the ones where bare Newton fails 18 of.
//! * `ported_step_rule_does_not_reopen_the_tail_crawl` -- the exponential-tail fixture stays
//!   `<= 40` iterations and strictly beats the strict-decrease variant (12 vs 119).
//!
//! `solve_v0` below is the pre-port loop, kept verbatim as the reference implementation those
//! guards measure against; sections A-F keep running the full comparison tables.
//!
//! ## What the old solve could not do (section F, 48 adversarial cases)
//!
//! On the ordinary grid the old solve converged everywhere: 84/84, residuals ~1e-16, agreeing
//! with the new root to 9.5e-13.  Its 11x cost came with no robustness dividend *visible on
//! that grid*.  Push the strike away from the money -- K = 100, or 1e6, against a bond worth
//! about 1 -- and the picture changes: **18 of 48 cases the old solve never converged**.  It
//! burned the full 50-iteration cap, and most of those ended on a `NaN` residual because a
//! Newton step overshot into exp overflow with no bracket to pull it back inside.  Worst old
//! root error: 3.5e1, i.e. 35,000bp from the reference.  The bracketed solve answered every
//! one of them.
//!
//! So the rewrite bought real robustness and paid 11x for it on the cases where the old solver
//! was already fine.  `v4` buys the robustness without the 11x.
//!
//! ## Layout
//!
//! | Module | What lives there |
//! |---|---|
//! | `fixtures` | curve shapes, coupon schedules, the strike ladder, the old-solver constants |
//! | `measurement` | the counted objective/derivative and the `Run` record every variant returns |
//! | `trace` | per-iteration trajectory reconstruction (`AB_TRACE=1`) |
//! | `variants` | the candidate solve loops `v2` / `v3` / `v4` and their shared bracket helper |
//! | `reference_loop` | the pre-port loop, kept as an independent reference implementation |
//! | `report` | medians, ratios, `Totals`, the old-vs-new verdict and the `compare` driver |
//! | `cost_sections` | sections A-D: per-case cost, aggregates, wall clock, the candidate fix |
//! | `exit_rules` | section E: exit-rule variants and the safety properties they must keep |
//! | `adversarial` | section F: the adversarial grid, where bare Newton breaks |
//! | `guards` | the always-on guards for the ported step rule |
//! | `threshold_scan` | the acceptance-ratio scan behind `MAX_MODEL_ERROR_RATIO` |
//!
//! Everything in here is test-only: the parent declares this module `#[cfg(test)]`, and the whole
//! module is additionally `#![cfg(test)]` so nothing here can reach a non-test build.

#![cfg(test)]

mod adversarial;
mod cost_sections;
mod exit_rules;
mod fixtures;
mod guards;
mod measurement;
mod reference_loop;
mod report;
mod threshold_scan;
mod trace;
mod variants;
