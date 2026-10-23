# geom.local_solve1d: Oscillation bracket search ignores failed samples (#481, 2026-10-23)

Unit: `geom.local_solve1d`; upstream defect WR-253; PR #668.

`LocalSolve1d::solve_oscillating` used the final attempted sample's F after its
search loop, even when that evaluation had failed. An evaluator that zeroes F
before returning an error therefore supplied a false sign change. This is the
same failure mode as the radius-of-curvature evaluator reported in #481.

The loop now records a sign change only after a successful evaluation. Failed
samples are still retried in the original order; successful bracket searches,
the best valid interior sample, existing one-sided brackets, the local solver's
accuracy and all formulas remain unchanged. `Unported` errors still propagate.

Two regressions failed before the correction: the existing cubic Newton-cycle
case with a zeroed final failed sample, and a search in which all eight interior
samples fail after writing zero or a value of the opposite sign. The corrected
cubic has the same no-root result, parameter, F and bracket flags as the clean
solve. A positive regression verifies that a valid later sample can still form
a bracket after an earlier evaluation fails, finding the root at 0.85.

The native test evaluator now supports zeroing F and F' on failure, to reproduce
the reported caller behavior independently in C++ and Rust. The fixed-size
14-value test ABI and existing failure modes 0 (leave outputs) and 1 (write F)
are preserved; mode 2 zeroes both outputs. The new saved fixture is
`oscillation_failed_last_guess_zeroed_outputs`; no solver reference code is
changed, no comparison removed and no tolerance widened.


The unchanged C++ solver is now confirmed natively on the zero-on-failure fixture:
`F=t^3-2t+2`, interval `[-3,3]`, initial guess `0`, and a failure exactly at
`0.94564800000197`. C++ returns `found=true`, `t=0.94564800000197`, F=0 and
`FoundAnswerConverged`; Rust returns `found=false`, `t=0.8316000000017324`,
actual F approximately `0.9119000984961294`, and `UnableToConverge`.
The failed sample's true F is about `0.95`. This comparison remains an explicit
known mismatch and the unit is partial, not green, until the reference is fixed.

The standalone macOS native run also reproduces the existing `gen1231` trig
mismatch (#654, already addressed by separate PR #655): C++ bracket F is
`3.493706235424554e-8` versus Rust `3.493706240975669e-8`. It is not caused by
the bracket-status fix and is not added to the known-mismatch list here.


Standalone Rust gates pass on macOS arm64 / Rust 1.98: 2,073 kernel tests,
formatting, workspace Clippy with and without `cpp` (`-D warnings`), wasm32
check, and default difftest tests. Native: 3,317 cases, 3,315 passes, one explicit
WR-253 known mismatch, no deferrals, and the one pre-existing #654 macOS trig
failure. The latter remains visible pending the separate PR #655 fix.


Combined compatibility validation (2026-10-24): main `d51978934` plus PRs
#642/#643/#644/#653/#655/#657/#660/#662/#668, unpublished verification commit
`ef437348d`: all 2,084 kernel tests, formatting, default/native workspace Clippy
with warnings denied and wasm32 check pass. All thirteen affected native units
have zero unexpected failures: five pass, eight remain partial with dependency
deferrals or explicit known mismatches. With the separate #655 trig correction,
`geom.local_solve1d` has 3,318 cases, 3,317 passes, one WR-253 known mismatch,
zero unexpected failures and zero deferrals. This is combined-source validation,
not a claim that the standalone macOS head or canonical Linux CI is green.
