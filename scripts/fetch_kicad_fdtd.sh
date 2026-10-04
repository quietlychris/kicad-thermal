#!/usr/bin/env bash
# Fetch the Python original, kicad_fdtd, at the commit this port follows, and apply the AMG preconditioner patch.
# Only needed to regenerate the test reference (scripts/gen_reference.py), to compare against Python
# (scripts/compare_live.py, bench/), or to export a board's geometry JSON (it needs KiCad's pcbnew).
#   scripts/fetch_kicad_fdtd.sh [DIR]      default DIR: third_party/kicad_fdtd (git-ignored)
# Then point the scripts at it with KICAD_FDTD=DIR (they default to third_party/kicad_fdtd).
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DIR=${1:-$ROOT/third_party/kicad_fdtd}
URL=https://github.com/Jerome-Graves/kicad_fdtd.git
COMMIT=74e21ebe6b87e19a7c709390d0a9661dc8574576
PATCH=$ROOT/patches/kicad_fdtd-amg-preconditioner.patch
if [ ! -d "$DIR/.git" ]; then
  mkdir -p "$(dirname "$DIR")"
  git clone --quiet "$URL" "$DIR"
fi
git -C "$DIR" fetch --quiet origin "$COMMIT" 2>/dev/null || true
git -C "$DIR" checkout --quiet --force "$COMMIT"
git -C "$DIR" reset --quiet --hard "$COMMIT"
git -C "$DIR" apply "$PATCH"
echo "kicad_fdtd $COMMIT with the AMG patch is at $DIR"
echo "Python dependencies: numpy scipy shapely pyamg (pyamg only for --precond amg/auto and gen_reference.py)"
