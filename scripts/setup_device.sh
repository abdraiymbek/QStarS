#!/bin/bash
# Device transfer (2026-10-29): first run on the new device after scripts/transfer_to_device.sh.
# usage: scripts/setup_device.sh [--check-only]
# Checks the toolchain and the copied runtime, rebuilds .venv from requirements-lock.txt (uv, Python 3.12),
# runs the test suite and reports the disk budget. --check-only: the checks alone (no venv, no tests).
set -uo pipefail
R=/Users/abraiymbek/Projects/QStarS
RT=/Users/abraiymbek/qsrd-runtime
cd "$R" || exit 1
export PATH="$HOME/.local/bin:/opt/homebrew/bin:$PATH"
fail=0
check() {   # check NAME COMMAND [HOW TO FIX]
  if eval "$2" >/dev/null 2>&1; then echo "ok    $1"; else echo "FAIL  $1${3:+ -> $3}"; fail=1; fi
}
check 'user abraiymbek, repo at the original path' '[ "$(whoami)" = abraiymbek ] && [ "$(pwd -P)" = "$R" ]'
check 'Apple Silicon' '[ "$(uname -m)" = arm64 ]'
check 'Xcode command line tools (linker, make, patch, git)' 'xcode-select -p' 'xcode-select --install'
check 'zstd (restoring archives)' 'command -v zstd' 'brew install zstd'
check 'uv (Python environment)' 'command -v uv' 'brew install uv'
check 'MESA SDK gfortran' "test -x $RT/mesasdk/bin/gfortran"
check 'MESA 24.08.1 tree' "test -d $RT/mesa/star && test -d $RT/mesa/data/eosDT_data"
for b in 17 22 31 32 33 34 35 36; do
  check "CMFGEN build $b tree" "test -d backends/cmfgen/build_$b/cur_cmf"
done
exe_libs() {   # every non-system library of the stage-1 executable exists (its @rpath is the MESA SDK)
  otool -L results/nk167_trap/exe_b36/cmfgen_dev.exe | awk 'NR>1{print $1}' | sed "s#@rpath#$RT/mesasdk#" |
    while read -r l; do case $l in /System/*|/usr/lib/*) ;; *) test -e "$l" || return 1 ;; esac; done
}
dangling() {   # symlinks of the working set (atomic data in sources/cache) that do not resolve
  find results sources -type l ! -exec test -e {} \; -print 2>/dev/null | head -1 | grep -q . && return 1 || return 0
}
check 'stage-1 executable (build 36) and its libraries' exe_libs
check 'stage-1 run directory and launch script' \
  'test -f results/nk167_trap/fixedT_b35/history.json && test -f results/nk167_trap/launch_fixedT_b36.sh'
check 'symlinks of the working set resolve' dangling
check 'archive index and records' 'test -f results/_archive/INDEX.json && ls results/_archive/record/*.tar.zst'
check 'git repository (remote reachable)' 'git status --short && git ls-remote origin HEAD' 'check the network; pushing needs gh auth login'
[ "${1:-}" = --check-only ] && { echo; du -sh "$R" "$RT"; exit $fail; }

echo "== Python environment (.venv from requirements-lock.txt)"
if [ ! -x .venv/bin/python ]; then
  uv venv --python 3.12 .venv || fail=1
fi
uv pip install --python .venv/bin/python -r requirements-lock.txt && \
  uv pip install --python .venv/bin/python --no-deps -e '.[test,audit]' || fail=1
echo "== test suite"
.venv/bin/python -m pytest -q 2>&1 | tail -3
[ "${PIPESTATUS[0]}" = 0 ] || fail=1
echo "== disk budget (target <= 100 GB for repo + runtime)"
du -sh "$R" "$RT"
[ $fail = 0 ] && echo "device ready" || echo "device NOT ready: fix the FAIL lines above and run again"
exit $fail
