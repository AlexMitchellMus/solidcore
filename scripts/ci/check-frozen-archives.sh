#!/usr/bin/env bash
# Fails if a retired shared archive, plan/findings.md or plan/upstream-issues.md, exists in
# HEAD's tree. Their records were split into one file each under plan/findings/ and
# plan/upstream-issues/ (scripts/split_archives.py) so parallel PRs never conflict
# (AGENTS.md section 5). Deleting them passes; re-adding either (e.g. a stale branch that
# appended to one) fails.
# Usage: check-frozen-archives.sh [<base sha>]   (the base is no longer needed; accepted so
# the CI call stays the same)
set -euo pipefail
fail=0
for f in plan/findings.md plan/upstream-issues.md; do
  if git cat-file -e "HEAD:$f" 2>/dev/null; then
    dir=${f%.md}
    echo "::error file=$f::$f was split into one file per record under $dir/ and must not come back. Put each new entry in its own file under $dir/ (see $dir/README.md); to update an old record, add an addendum file there."
    fail=1
  fi
done
exit $fail
