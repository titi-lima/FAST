#!/usr/bin/env bash
set -euo pipefail

# One-shot wrapper: build, generate bbox, prioritize, map to selectors, run in order

ROOT_DIR="$(cd "$(dirname "$0")"/../.. && pwd)"
EX_DIR="$ROOT_DIR/ant"    # respect user's possible path changes
DATA_DIR="$EX_DIR/.fast/in"
OUT_DIR="$EX_DIR/.fast/out"
BUILD_XML="$EX_DIR/build.xml"
LIB_JAR="$EX_DIR/lib/junit-platform-console-standalone.jar"

ALG="FAST-pw"
REPS=3

mkdir -p "$DATA_DIR" "$OUT_DIR"

echo "[1/5] Ensure JUnit console present"
if [ ! -f "$LIB_JAR" ]; then
  "$EX_DIR/scripts/get-junit.sh"
fi

echo "[2/5] Build (Java 21)"
ant -f "$BUILD_XML" clean compile | cat

echo "[3/5] Generate bbox from tests (selectors + tokens)"
python3 "$EX_DIR/scripts/tests-map.py" \
  --tests-dir "$EX_DIR/src/test/java" \
  --out "$DATA_DIR/selectors.txt"

python3 "$EX_DIR/scripts/generate-bbox.py" \
  --tests-dir "$EX_DIR/src/test/java" \
  --out "$DATA_DIR/ant-sample-bbox.txt"

echo "[4/5] Prioritize with FAST TCP ($ALG, reps=$REPS)"
python3 -m fast_tcp.cli \
  --test-dir "$DATA_DIR" \
  --algo "$ALG" \
  --entity bbox \
  --repetitions "$REPS" \
  --file-naming entity-suffix \
  --output-dir "$OUT_DIR" | cat

echo "[5/5] Map prioritized IDs to selectors and run tests in that order"
python3 "$EX_DIR/scripts/build-prioritized-selectors.py" \
  --out-dir "$OUT_DIR" \
  --selectors "$DATA_DIR/selectors.txt" \
  --out "$DATA_DIR/prioritized-selectors.txt"

# Build classpath
MAIN_CP="$EX_DIR/build/classes/main"
TEST_CP="$EX_DIR/build/classes/test"
CP="$MAIN_CP:$TEST_CP:$LIB_JAR"

ARGS=()
while IFS= read -r line; do
  [ -n "$line" ] && ARGS+=("--select-method" "$line")
done < "$DATA_DIR/prioritized-selectors.txt"

echo "Running JUnit Console in prioritized order..."
java -cp "$CP" org.junit.platform.console.ConsoleLauncher \
  --fail-if-no-tests \
  --class-path "$MAIN_CP:$TEST_CP" \
  "${ARGS[@]}"

echo "Done."


