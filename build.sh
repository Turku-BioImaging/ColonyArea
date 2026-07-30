#!/usr/bin/env bash
set -euo pipefail

# Note: set the IJ_VAR and VERSION environment variables before running this script.

# Path to ij-1.54p.jar from your local Fiji installation
IJ_JAR="${IJ_JAR:-/path/to/Fiji/jars/ij-1.54p.jar}"
# Plugin version appended to the JAR and zip filenames (e.g. IJ_JAR=... VERSION=1.1.0 ./build.sh)
VERSION="${VERSION:-0.0.0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
OUT_DIR="$BUILD_DIR/classes"
DIST_DIR="$BUILD_DIR/dist"

echo "Cleaning build directory..."
rm -rf "$BUILD_DIR"
mkdir -p "$OUT_DIR" "$DIST_DIR/ColonyArea"

echo "Compiling ColonyArea_.java..."
javac --release 21 \
  -cp "$IJ_JAR" \
  -d "$OUT_DIR" \
  "$SCRIPT_DIR/src/ColonyArea_.java"

echo "Packaging ColonyArea_-${VERSION}.jar..."
jar cf "$DIST_DIR/ColonyArea_-${VERSION}.jar" -C "$OUT_DIR" .

echo "Copying macro files..."
cp "$SCRIPT_DIR/src/Colony_measurer.ijm"        "$DIST_DIR/ColonyArea/"
cp "$SCRIPT_DIR/src/Colony_thresolder.ijm"       "$DIST_DIR/ColonyArea/"
cp "$SCRIPT_DIR/src/Manual_colony_thresholder.ijm" "$DIST_DIR/ColonyArea/"

echo "Creating zip archive..."
cd "$DIST_DIR"
zip -r "$BUILD_DIR/ColonyArea-${VERSION}.zip" "ColonyArea_-${VERSION}.jar" ColonyArea/

echo "Done. Output: $BUILD_DIR/ColonyArea-${VERSION}.zip"