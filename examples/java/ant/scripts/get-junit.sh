#!/usr/bin/env bash
set -euo pipefail

# Download JUnit Platform Console Standalone to lib/
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ANT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
LIB_DIR="$ANT_DIR/lib"
mkdir -p "$LIB_DIR"

VERSION="1.10.3"
JAR="junit-platform-console-standalone-${VERSION}.jar"
URL="https://repo1.maven.org/maven2/org/junit/platform/junit-platform-console-standalone/${VERSION}/${JAR}"

echo "Downloading $JAR ..."
curl -fsSL "$URL" -o "$LIB_DIR/$JAR"
ln -sf "$JAR" "$LIB_DIR/junit-platform-console-standalone.jar"
echo "Saved to $LIB_DIR/$JAR"


