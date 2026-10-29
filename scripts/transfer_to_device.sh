#!/bin/bash
# Device transfer (2026-10-29): copy the QStarS working set to an always-on Apple Silicon Mac with the same user
# (abraiymbek) and the same absolute paths, then run scripts/setup_device.sh there.
#
# usage: scripts/transfer_to_device.sh HOST [--dry-run] [--verify]
#   HOST      ssh destination, e.g. abraiymbek@mac-mini.local (Remote Login on; this Mac's key authorized there)
#   --dry-run list and size what would be copied, copy nothing
#   --verify  after copying, compare every file by checksum (reads both sides once more)
#
# Copied to the same paths: the repo with .git; results/ (the stage-1 run, its executables and inputs, and
# results/_archive record + manifests + INDEX); CMFGEN builds 17, 22, 31-36 (tests and the current executable);
# MESA/FreeEOS builds; sources/cache without the two archives whose extracted trees are copied; ~/qsrd-runtime mesa,
# mesasdk; Claude's project notes. Stay on this Mac only: results/_archive/full (the lossless archives), other CMFGEN
# builds, the archives sources/cache/atomic_data_21jun23.tar.gz and three_recent_models.tar.gz, the .venv (rebuilt
# there from requirements-lock.txt) and the old MESA run directories in ~/qsrd-runtime.
set -euo pipefail
HOST=${1:?usage: transfer_to_device.sh HOST [--dry-run] [--verify]}
shift
DRY=; VERIFY=
for a in "$@"; do
  case $a in --dry-run) DRY=-n ;; --verify) VERIFY=1 ;; *) echo "unknown option $a" >&2; exit 2 ;; esac
done
R=/Users/abraiymbek/Projects/QStarS
RT=/Users/abraiymbek/qsrd-runtime
CL=/Users/abraiymbek/.claude/projects/-Users-abraiymbek-Projects-QStarS
cd "$R"

raw=$(.venv/bin/python scripts/archive_results.py plan | grep -c '(raw)' || true)
if [ "$raw" != 0 ]; then
  echo "results/ still holds $raw unpacked units (scripts/archive_results.py plan): pack them first" >&2; exit 1
fi

if [ "$HOST" = local ]; then   # --dry-run sizing without a device: destination is a fresh scratch directory
  DEST=${TMPDIR:-/tmp}/qstars_transfer_dryrun; RSH=; PRE=
  [ -n "$DRY" ] || { echo "local is for --dry-run only" >&2; exit 2; }
else
  info=$(ssh -o BatchMode=yes -o ConnectTimeout=10 "$HOST" \
    'echo "$(whoami) $(uname -m) $(sw_vers -productVersion) $(df -g ~ | awk "NR==2{print \$4}")"; test -d ~/Projects && echo projects')
  set -- $info
  [ "$1" = abraiymbek ] || { echo "remote user is $1, not abraiymbek: the absolute paths would break" >&2; exit 1; }
  [ "$2" = arm64 ] || { echo "remote is $2, not Apple Silicon: the native builds would not run" >&2; exit 1; }
  [ "${5:-}" = projects ] || { echo "remote ~/Projects is missing" >&2; exit 1; }
  echo "remote: $1, $2, macOS $3, $4 GB free"
  [ "$4" -ge 60 ] || { echo "need at least 60 GB free on the remote" >&2; exit 1; }
  ssh "$HOST" "mkdir -p $R $RT $CL"
  DEST=; PRE="$HOST:"
fi

sync() {   # sync SRC DEST_DIR [rsync filter args...]
  local src=$1 dst=$2; shift 2
  [ -z "$PRE" ] && mkdir -p "$DEST$dst"
  rsync -a $DRY --partial --stats "$@" "$src" "$PRE$DEST$dst" | grep -E 'Number of files|Total file size|Total transferred|^sent' || true
}

echo "== repository and working set"
sync "$R/" "$R/" \
  --exclude=/.venv/ --exclude=/.pytest_cache/ --exclude=__pycache__/ --exclude=.DS_Store \
  --exclude=/results/_archive/full/ --exclude='/results/_archive/*.part' --exclude='/results/_archive/*/*.part' \
  --include=/backends/cmfgen/build_17/ --include=/backends/cmfgen/build_22/ \
  --include='/backends/cmfgen/build_3[1-6]/' --exclude='/backends/cmfgen/build_*/' \
  --exclude=/sources/cache/atomic_data_21jun23.tar.gz --exclude=/sources/cache/three_recent_models.tar.gz
echo "== MESA and its SDK"
sync "$RT/mesa" "$RT/"
sync "$RT/mesasdk" "$RT/"
echo "== Claude project notes"
sync "$CL/" "$CL/"

[ -n "$DRY" ] && { echo "dry run: nothing copied"; exit 0; }
if [ -n "$VERIFY" ]; then
  echo "== checksum verification"
  diff=$( { rsync -anc --out-format='%n' --exclude=/.venv/ --exclude=/.pytest_cache/ --exclude=__pycache__/ \
      --exclude=.DS_Store --exclude=/results/_archive/full/ --exclude='/results/_archive/*.part' \
      --exclude='/results/_archive/*/*.part' --include=/backends/cmfgen/build_17/ \
      --include=/backends/cmfgen/build_22/ --include='/backends/cmfgen/build_3[1-6]/' \
      --exclude='/backends/cmfgen/build_*/' --exclude=/sources/cache/atomic_data_21jun23.tar.gz \
      --exclude=/sources/cache/three_recent_models.tar.gz "$R/" "$HOST:$R/"
    rsync -anc --out-format='%n' "$RT/mesa" "$RT/mesasdk" "$HOST:$RT/"; } | grep -v '/$' || true)
  if [ -n "$diff" ]; then
    echo "$diff" | head -20; echo "$(echo "$diff" | wc -l) files differ: run the transfer again" >&2; exit 1
  fi
  echo "every copied file identical by checksum"
fi
echo "== setup on the device"
ssh "$HOST" "cd $R && /bin/bash scripts/setup_device.sh"
