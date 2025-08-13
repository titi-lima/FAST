#!/usr/bin/env bash
set -euo pipefail

# End-to-end: generate bbox, prioritize, map to test names, build, run prioritized order

ROOT_DIR="$(cd "$(dirname "$0")"/../.. && pwd)"   # examples/go
EX_DIR="$ROOT_DIR/gotest"                            # this example root
REPO_DIR="$(cd "$ROOT_DIR"/../.. && pwd)"          # repo root
DATA_DIR="$EX_DIR/.fast/in"
OUT_DIR="$EX_DIR/.fast/out"

ALG="FAST-pw"
REPS=3

mkdir -p "$DATA_DIR" "$OUT_DIR"

PY="$REPO_DIR/virtual/bin/python"
if [ ! -x "$PY" ]; then
  PY="python3"
fi

echo "[1/5] Generate bbox (tokens) and test names"
"$PY" "$EX_DIR/scripts/tests-map.py" \
  --tests-dir "$EX_DIR" \
  --out "$DATA_DIR/test-names.txt"

"$PY" "$EX_DIR/scripts/generate-bbox.py" \
  --tests-dir "$EX_DIR" \
  --out "$DATA_DIR/gotest-bbox.txt"

echo "[2/5] Prioritize with FAST TCP ($ALG, reps=$REPS)"
"$PY" -m fast_tcp.cli \
  --test-dir "$DATA_DIR" \
  --algo "$ALG" \
  --entity bbox \
  --repetitions "$REPS" \
  --file-naming entity-suffix \
  --output-dir "$OUT_DIR" | cat

echo "[3/5] Map prioritized IDs to Go test names"
"$PY" "$EX_DIR/scripts/build-prioritized-tests.py" \
  --out-dir "$OUT_DIR" \
  --tests "$DATA_DIR/test-names.txt" \
  --out "$DATA_DIR/prioritized-tests.txt"

echo "[4/5] Build test binary"
cd "$EX_DIR"

BIN_DIR="$EX_DIR/.fast/bin"
mkdir -p "$BIN_DIR"

# Build a test binary per package directory
while IFS= read -r d; do
  [ -z "$d" ] && continue
  (cd "$d" && pkg_name=$(basename "$d") && go test -c -o "$BIN_DIR/${pkg_name}.test" || true)
done < <(go list -f '{{.Dir}}' ./...)

echo "[5/5] Run tests in prioritized order"
while IFS= read -r test_name; do
  [ -z "$test_name" ] && continue
  # Run across binaries until one matches the test
  matched=0
  for bin in "$BIN_DIR"/*.test; do
    if "$bin" -test.list . | grep -E "^$test_name$" >/dev/null 2>&1; then
      echo "==> $test_name"
      "$bin" -test.run "^$test_name$" | cat || true
      matched=1
      break
    fi
  done
  if [ "$matched" -eq 0 ]; then
    echo "Warning: test $test_name not found in built binaries"
  fi
done < "$DATA_DIR/prioritized-tests.txt"

echo "Done."


