#!/usr/bin/env bash
set -eu

SRC_DIR="${SRC_DIR:-$(pwd)}"
BUILD_ROOT="${BUILD_ROOT:-$SRC_DIR/_build_matrix_cmake}"
LOG_DIR="$BUILD_ROOT/logs"
CTEST_CMD="${CTEST:-ctest}"

mkdir -p "$LOG_DIR"

exec "$CTEST_CMD" \
  -S "$SRC_DIR/qa/cmake_matrix.cmake" \
  -DSRC_DIR="$SRC_DIR" \
  -DBUILD_ROOT="$BUILD_ROOT" \
  -DLOG_DIR="$LOG_DIR"
