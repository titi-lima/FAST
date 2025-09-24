#!/usr/bin/env bash
set -euo pipefail

# End-to-end: generate bbox, prioritize, map to test names, build, run prioritized order

PROJECT_DIR="$(pwd)"
DATA_DIR="$PROJECT_DIR/.fast/in"
OUT_DIR="$PROJECT_DIR/.fast/out"
BIN_DIR="$PROJECT_DIR/.fast/bin"

ALG="FAST-pw"
REPS=3

mkdir -p "$DATA_DIR" "$OUT_DIR" "$BIN_DIR"

PY="python3"

echo "[1/5] Generate bbox (tokens) and test names"
"$PY" "$PROJECT_DIR/.fast/tools/gotest/tests-map.py" \
  --tests-dir "$PROJECT_DIR" \
  --out "$DATA_DIR/test-names.txt"

"$PY" "$PROJECT_DIR/.fast/tools/gotest/generate-bbox.py" \
  --tests-dir "$PROJECT_DIR" \
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
"$PY" "$PROJECT_DIR/.fast/tools/gotest/build-prioritized-tests.py" \
  --out-dir "$OUT_DIR" \
  --tests "$DATA_DIR/test-names.txt" \
  --out "$DATA_DIR/prioritized-tests.txt"

echo "[4/5] Build test binaries"
while IFS= read -r d; do
  [ -z "$d" ] && continue
  (cd "$d" && pkg_name=$(basename "$d") && go test -c -o "$BIN_DIR/${pkg_name}.test" || true)
done < <(go list -f '{{.Dir}}' ./...)

echo "[5/5] Run tests in prioritized order"
while IFS= read -r test_name; do
  [ -z "$test_name" ] && continue
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


