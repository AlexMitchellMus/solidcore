Updates: WR-253

#### Issue #481 / WR-253 follow-up. Failed oscillation sample confirmed as a false root
- **Severity:** Wrong result
- **Area:** Topology, attributes, solvers

- **Location:** `source/SMLib/src/SmLocalSolve1d.cpp:470-495`, the sign test after the interior search.
- **Trigger:** `F=t^3-2t+2` on `[-3,3]`, guess zero, with an evaluator that zeroes its outputs and fails at the final interior sample `0.94564800000197`.
- **Confirmed:** yes, now reproduced through the unchanged native solver in the differential harness. C++ reports that failed sample as a converged root with F=0, although its actual F is about 0.95; the corrected Rust search reports no root.
- **Rust port:** PR #668 requires a successfully evaluated sign change before installing an interior bracket. Failed samples still allow later valid samples and fallback search. Three direct regressions cover these paths; the native defect is one explicit known mismatch, keeping the unit partial. This supersedes WR-253's earlier no-native-confirmation/faithful-port status without rewriting its historical entry.
