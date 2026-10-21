#!/usr/bin/env bash
# Fetches the difftest reports (reports/difftest/<unit>.json) of main's latest CI run.
#
# Reports are not committed: the CI difftest job on every push to main uploads them as
# the artifact `difftest-reports`. This downloads the newest one into reports/difftest/
# (replacing its *.json). Needs `gh` with read access to Actions (GH_TOKEN in CI).
# When no artifact is available (none uploaded yet, or expired) it falls back to the
# last committed set in git history (needs a non-shallow clone), and exits 1 if that
# fails too. Running difftest locally regenerates the reports of the units you run.
#
#   scripts/fetch-reports.sh [DEST]          # DEST defaults to reports/difftest
set -euo pipefail

root=$(git rev-parse --show-toplevel)
dest=${1:-$root/reports/difftest}
repo=${GITHUB_REPOSITORY:-tangency-labs/solidcore}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

install() {
  local src=$1 n
  n=$(find "$src" -maxdepth 1 -name '*.json' | wc -l)
  [ "$n" -gt 0 ] || return 1
  mkdir -p "$dest"
  find "$dest" -maxdepth 1 -name '*.json' -delete
  find "$src" -maxdepth 1 -name '*.json' -exec mv -t "$dest" {} +
  echo "fetch-reports: $n reports -> ${dest#"$root"/} ($2)"
}

if command -v gh >/dev/null; then
  # Newest unexpired artifacts first; only ones uploaded by a run on main.
  ids=$(gh api "repos/$repo/actions/artifacts?name=difftest-reports&per_page=20" \
          -q '.artifacts[] | select(.expired | not) | select(.workflow_run.head_branch == "main")
              | .workflow_run.id' 2>/dev/null | head -n 5 || true)
  for id in $ids; do
    if gh run download "$id" -R "$repo" -n difftest-reports -D "$tmp/run" >/dev/null 2>&1 \
        && install "$tmp/run" "CI run $id"; then
      exit 0
    fi
    rm -rf "$tmp/run"
  done
  echo "fetch-reports: no difftest-reports artifact from main's CI" >&2
fi

# Fallback: the set committed before reports moved to CI artifacts.
last=$(git -C "$root" rev-list -1 HEAD -- reports/difftest 2>/dev/null || true)
if [ -n "$last" ]; then
  for c in "$last" "$last^"; do
    if git -C "$root" cat-file -e "$c:reports/difftest" 2>/dev/null; then
      mkdir -p "$tmp/git"
      git -C "$root" archive "$c" reports/difftest | tar -x -C "$tmp/git"
      install "$tmp/git/reports/difftest" "committed in ${c}" && exit 0
    fi
  done
fi
echo "fetch-reports: no reports found" >&2
exit 1
