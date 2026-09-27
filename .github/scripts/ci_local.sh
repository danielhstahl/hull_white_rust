#!/usr/bin/env bash
#
# Run the CI jobs locally, in the order the workflows run them, and print one
# PASS/FAIL/SKIP line per job.
#
# This is the branch-level dry run that stands in for a GitHub-hosted run when you
# cannot (or do not want to) push: the same commands each job runs, with the same
# exit codes, on this machine.  It does not replace CI; it answers "would these
# jobs be green, and is each one's status meaningful?" without a push.
#
# Usage:
#   bash .github/scripts/ci_local.sh                 # every job
#   bash .github/scripts/ci_local.sh fmt clippy      # a named subset
#   bash .github/scripts/ci_local.sh --list          # job names
#   bash .github/scripts/ci_local.sh --bench         # actually run the benches
#                                                  # (compile-only by default)
#
# Nightly jobs need a nightly toolchain: `rustup toolchain install nightly --component clippy`.
# Coverage needs `cargo-llvm-cov`.  A missing tool is reported SKIP, not PASS --
# a job that never ran must never look like a job that passed.

# No `-e`: the point is to run everything and report, not to die at the first red.
set -uo pipefail

# Every command below is relative to the crate root; make that true regardless of
# where the script was invoked from.
cd "$(git rev-parse --show-toplevel)" || exit 1

run_benches=0
jobs=()
for arg in "$@"; do
  case "$arg" in
    --list)
      echo "fmt clippy clippy-all-targets build test test-slow coverage doc bench publish"
      exit 0
      ;;
    --bench) run_benches=1 ;;
    *) jobs+=("$arg") ;;
  esac
done
if [ "${#jobs[@]}" -eq 0 ]; then
  jobs=(fmt clippy clippy-all-targets build test test-slow coverage doc bench publish)
fi

names=()
results=()
notes=()
overall=0

record() { # name, status, note
  names+=("$1")
  results+=("$2")
  notes+=("$3")
  if [ "$2" = "FAIL" ]; then overall=1; fi
}

run() { # name, required-command-or-empty, command...
  local name="$1"
  shift
  local requires="$1"
  shift
  if [ -n "$requires" ] && ! command -v "$requires" >/dev/null 2>&1; then
    echo "-- $name: SKIP (needs '$requires', not installed)"
    record "$name" "SKIP" "needs $requires"
    return 0
  fi
  echo "-- $name: $*"
  local start now code
  start=$(date +%s)
  if "$@"; then
    now=$(( $(date +%s) - start ))
    echo "-- $name: PASS (${now}s)"
    record "$name" "PASS" "${now}s: $*"
  else
    code=$?
    now=$(( $(date +%s) - start ))
    echo "-- $name: FAIL (exit ${code}, ${now}s)"
    record "$name" "FAIL" "exit ${code}: $*"
  fi
}

nightly_ready() {
  command -v rustup >/dev/null 2>&1 && rustup run nightly cargo --version >/dev/null 2>&1
}

for job in "${jobs[@]}"; do
  case "$job" in
    # lint.yml / fmt
    fmt)
      run "fmt" cargo cargo fmt --all -- --check
      ;;
    # lint.yml / clippy (stable).  Not `--all-targets`: that pulls in benches/,
    # which need nightly's `#![feature(test)]`.
    clippy)
      run "clippy" cargo cargo clippy --lib --bins --tests --all-features -- -D warnings
      ;;
    # lint.yml / clippy --all-targets (nightly).  Adds the bench target.
    clippy-all-targets)
      if nightly_ready; then
        run "clippy-all-targets" "" cargo +nightly clippy --all-targets --all-features -- -D warnings
      else
        echo "-- clippy-all-targets: SKIP (no nightly toolchain: rustup toolchain install nightly)"
        record "clippy-all-targets" "SKIP" "no nightly toolchain"
      fi
      ;;
    build)
      run "build" cargo cargo build --verbose
      ;;
    test)
      run "test" cargo cargo test --verbose
      ;;
    test-slow)
      run "test-slow" cargo cargo test --verbose --features slow
      ;;
    coverage)
      if command -v cargo-llvm-cov >/dev/null 2>&1; then
        run "coverage" "" cargo llvm-cov --all-features --workspace --summary-only
      else
        echo "-- coverage: SKIP (cargo-llvm-cov not installed)"
        record "coverage" "SKIP" "cargo-llvm-cov not installed"
      fi
      ;;
    doc)
      run "doc" cargo cargo doc --no-deps
      ;;
    bench)
      # Compile-only by default: the benches take minutes and their timings are
      # meaningless on shared/virtualised hardware.  `--bench` runs them anyway.
      if nightly_ready; then
        if [ "$run_benches" -eq 1 ]; then
          run "bench" "" cargo +nightly bench --features test-support
        else
          run "bench" "" cargo +nightly bench --features test-support --no-run
        fi
      else
        echo "-- bench: SKIP (no nightly toolchain)"
        record "bench" "SKIP" "no nightly toolchain"
      fi
      ;;
    publish)
      # The dry run needs the tracked tree clean, because that is what the real job
      # gets: a fresh checkout.  Uncommitted local work is not a red signal about
      # the branch, so it reports SKIP with the reason instead of a failure.
      if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
        echo "-- publish: SKIP (tracked tree is dirty; the real job runs on a clean checkout)"
        record "publish" "SKIP" "tracked tree dirty"
      else
        run "publish" bash bash .github/scripts/publish.sh --dry-run
      fi
      ;;
    *)
      echo "-- $job: unknown job (see --list)" >&2
      record "$job" "FAIL" "unknown job"
      ;;
  esac
done

echo
echo "=================================================================="
printf '%-20s %-6s %s\n' "JOB" "STATUS" "DETAIL"
echo "------------------------------------------------------------------"
i=0
while [ "$i" -lt "${#names[@]}" ]; do
  printf '%-20s %-6s %s\n' "${names[$i]}" "${results[$i]}" "${notes[$i]}"
  i=$((i + 1))
done
echo "=================================================================="
if [ "$overall" -eq 0 ]; then
  echo "Overall: every job that ran is green"
else
  echo "Overall: at least one job FAILED"
fi
exit "$overall"
