#!/bin/bash
# QStarS lives on this T9 drive (moved 2026-10-30): two APFS sparse-bundle images, each attached at the project's
# original path on the Mac mini, so every absolute path in the code, run records and executables stays valid.
#   QStarS.sparsebundle        -> /Users/abraiymbek/Projects/QStarS   (repository, .venv, results/)
#   qsrd-runtime.sparsebundle  -> /Users/abraiymbek/qsrd-runtime      (MESA, MESA SDK; the executables' @rpath)
# Run after every restart of the Mac (and after reconnecting the T9) before any QStarS work. Idempotent.
# usage: /Volumes/T9/QStarS/attach_qstars.sh
HERE=$(cd "$(dirname "$0")" && pwd)
status=0
attach() {   # image mountpoint marker
  if [ -f "$2/$3" ]; then
    echo "attached: $2"
    return 0
  fi
  [ -d "$HERE/$1" ] || { echo "image not found: $HERE/$1"; return 1; }
  if [ -e "$2" ] && [ -n "$(ls -A "$2" | grep -v '^README_NOT_ATTACHED.txt$')" ]; then
    echo "refusing to attach over $2: it is not an empty mount point"
    return 1
  fi
  mkdir -p "$2"
  hdiutil attach -nobrowse -mountpoint "$2" "$HERE/$1" > /dev/null || return 1
  if [ -f "$2/$3" ]; then echo "attached: $2"; else echo "attached $1 at $2, but marker $3 is missing"; return 1; fi
}
attach QStarS.sparsebundle /Users/abraiymbek/Projects/QStarS .qstars_project_volume || status=1
attach qsrd-runtime.sparsebundle /Users/abraiymbek/qsrd-runtime .qstars_runtime_volume || status=1
exit $status
