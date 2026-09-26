# Changelog

## 0.9.0 — one curve object (breaking)

### The decision

`HullWhite<'a, T, U>` took **two** closures — a cumulative `yield_curve` and an
`forward_curve` — that are one thing. The forward is derivable from the yield
(`f(0,t) = d/dt y(t)` under this crate's convention), so nothing stopped a caller
from passing a pair that disagreed, and a model built on a mutually inconsistent
pair prices every instrument wrong with no diagnostic. The generic pair also leaked
into user code: every function that took a model had to name two curve closure
types it did not care about.

The chosen shape is **a trait, `YieldCurve`, with the cumulative yield as its one
required method and the forward derived from it** — plus a construction-time
consistency check for anybody who supplies the forward explicitly:

```rust
pub trait YieldCurve: Sync {
    fn zero_yield(&self, t: f64) -> f64;                     // required: y(t)
    fn discount(&self, t: f64) -> f64 { (-self.zero_yield(t)).exp() }   // P(0,t)
    fn forward(&self, t: f64) -> f64 { /* 5-point central difference */ }
}
```

Why this shape rather than the alternatives:

* **One required method** because the cumulative yield is the actual primitive:
  `P(0,t) = exp(-y(t))` and `f(0,t) = y'(t)`. A caller with only a yield curve
  gets a working model with nothing else to write, and there is then no second
  number for it to disagree with.
* **`forward` overridable** because a finite difference is not always wanted.
  Every closed-form curve in this crate — the test fixtures, a Vasicek calibration,
  a bootstrapped piecewise curve — knows its own derivative, and overriding costs
  one method and takes the finite difference out of the pricing path.
* **A consistency check anyway**, rather than "derive it and there is nothing to
  check", because the case that matters is the one where two numbers are supplied:
  a caller migrating from `init`, or a curve implemented from a calibration where
  the yield and the forward were fitted separately. The check cannot help someone
  who derived the forward; it exists precisely for the people who did not.
* **Dynamic dispatch over generics.** `HullWhite` holds a `&dyn YieldCurve`, so the
  type parameter goes away entirely. The curve evaluations are already two indirections
  deep inside a Jamshidian solve that spends ~1,400 evaluations per option; the
  loss against a monomorphised closure is not measurable against that, and the type
  a consumer writes is `HullWhite<'a>`.

### The convention, stated

`zero_yield` is **cumulative**: `y(t) = integral of f(0,.) over [0, t]`, so
`P(0,t) = exp(-y(t))` with no division by maturity. The consistency relation under
that convention is the plain `f = y'`.

This is worth spelling out because the *annualised* zero rate `z(t) = y(t)/t` — what
a quoting screen prints — satisfies a different relation to the same forward,
`f(t) = z(t) + t z'(t)`. A `YieldCurve` implemented against the annualised zero
is not "off by a factor of t" in some distant price: it fails the construction
check, by name, at the first probe.

### The check

`HullWhite::new` compares `curve.forward(t)` against the 5-point finite-difference
derivative of `curve.zero_yield(t)` at every time in `CURVE_PROBE_TIMES`
(`0.25, 0.5, 1.0, 2.0, 5.0, 10.0` years), and rejects the curve if any probe
disagrees by more than

```
tolerance = max(1e-6, 1e-4 * |forward(t)|)
```

i.e. a hundredth of a basis point, relaxing to one part in 10,000 of the forward
where the forward is large. The failure is
`HullWhiteError::InvalidInput` carrying the worst probe time, the supplied forward,
the derived one, the error and the tolerance. A curve that cannot be evaluated at
any probe time is rejected too — "nothing was checked" is not a pass.

A supplied forward that is **not finite** over a finite yield is always rejected, and
the relative term is deliberately not applied to it: `1e-4 * |inf|` is `inf`, and
`inf > inf` is false, so growing the tolerance by the magnitude of the supplied
forward would certify `forward(t) = infinity` as consistent and then carry that
infinity into `phi(t)`, `mu_r` and the Jamshidian bracket. For a non-finite
`forward` the tolerance stays at the `1e-6` absolute floor, so the comparison is
`inf > 1e-6` — true — and the curve is refused. `NaN` was already caught by this
rule; `±infinity` was the case that slipped, and it is the nastier one because a
`NaN` price announces itself and an infinite drift does not.

`0.0` is deliberately not a probe time: the derivative there is one-sided and the
least accurate, and the front of the curve already has its own gate — a curve whose
forward is not finite at `0` cannot supply `r(0)` and is reported by
`HullWhite::short_rate_now`, which is a different fact about a different thing.

Measured against the crate's own fixtures the check has six orders of headroom
(worst `2.0e-12` on `high_vol`, tolerance `5.1e-6`), so it fires on curve-shape
errors and not on rounding.

### Breaking changes

| Before 0.9 | 0.9 |
|---|---|
| `HullWhite<'a, T, U>` | `HullWhite<'a>` |
| `HullWhite::init(a, sigma, &yield_curve, &forward_curve)` | `HullWhite::new(a, sigma, &curve)` |
| `fn f<T: Fn(f64)->f64 + Sync, U: Fn(f64)->f64 + Sync>(m: &HullWhite<'_, T, U>)` | `fn f(m: &HullWhite)` |
| two closures, trusted | one curve, checked |

`HullWhite::init` is `#[deprecated(since = "0.9.0")]`. It stays compiling and
deprecated for the whole of `0.9.x`, and is **removed in `0.10.0`** — the breaking
type collapse ships with `0.9.0`, the removal trails it by one release so a caller
can adopt the new constructor without the two changing underneath at once.
**The deprecated call is not behaviour-preserving**: it wraps its two closures in
`from_yield_and_forward` and runs the same consistency check, so a pair that used
to price (wrongly) is now an error. Migrate off it rather than silencing the
deprecation.

### Migration

One closure — the cumulative yield. The forward is derived:

```rust
// 0.8
let y = |t: f64| 0.05 * t + 0.01 * t * t;
let f = |t: f64| 0.05 + 0.02 * t;
let m = HullWhite::init(a, sigma, &y, &f)?;

// 0.9 -- you already have everything the model needs in `y`
let curve = hull_white::from_yield(|t: f64| 0.05 * t + 0.01 * t * t);
let m = HullWhite::new(a, sigma, &curve)?;
```

Both halves in closed form — the literal migration, no finite difference:

```rust
// 0.9
let curve = hull_white::from_yield_and_forward(y, f);
let m = HullWhite::new(a, sigma, &curve)?;
```

Your own curve — a bootstrapped discount curve, a piecewise-linear forward, a
calibration object. Implement the trait; override `forward` if the derivative is
known and `discount` if the curve is stored as prices rather than yields:

```rust
struct MyCurve { /* ... */ }
impl YieldCurve for MyCurve {
    fn zero_yield(&self, t: f64) -> f64 { /* ... */ }
    fn forward(&self, t: f64) -> f64 { /* analytic derivative */ }
}
```

`HullWhite::curve()` returns the curve back (`&dyn YieldCurve`), so
`model.curve().discount(t)` is available wherever a price came from it.

### If the migration check rejects a pair that used to work

It is telling you the pair was inconsistent, which is to say the prices were
already wrong and nothing had said so. Three ways out, in order of preference:

1. Supply the actual derivative of the yield you passed. This is the fix.
2. Drop the forward and use `from_yield(y)` — the model derives it, and there is
   then no pair to be inconsistent.
3. Implement `YieldCurve` on your own type and make the two accessors agree; the
   error message names the probe time and both values, so the divergence is
   findable.

There is deliberately no "skip the check" flag. A model that will not build on an
inconsistent curve is the feature.

### Also in 0.9.0

* `pub mod curves` is public: `YieldCurve`, `from_yield`, `from_yield_and_forward`,
  `ForwardInconsistency`, `max_forward_inconsistency`, `validate_curve`,
  `CURVE_PROBE_TIMES`, `FORWARD_CONSISTENCY_TOLERANCE`,
  `FORWARD_CONSISTENCY_RELATIVE_TOLERANCE`, `forward_consistency_tolerance`.
* `HullWhite` implements `Debug` (model parameters; the curve is a trait object
  and is named rather than dumped).
* `test_support::hw_curves` / `Scenario::curves`, which returned a
  `(yield_closure, forward_closure)` tuple, are replaced by `HwCurve` /
  `hw_curve` / `Scenario::curve` returning the single object. Test-only and
  `#[doc(hidden)]`, but benches building on `--features test-support` need the
  rename.
* The test-suite fixtures are unchanged in value: every price the suite asserts is
  the same number it was before the refactor, rewired through the new constructor.
