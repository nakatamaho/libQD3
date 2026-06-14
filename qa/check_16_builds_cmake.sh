#!/usr/bin/env bash
set -eu

SRC_DIR="${SRC_DIR:-$(pwd)}"
BUILD_ROOT="${BUILD_ROOT:-$SRC_DIR/_build_matrix_cmake}"
LOG_DIR="$BUILD_ROOT/logs"
CTEST_CMD="${CTEST:-ctest}"

# Set ENABLE_MPFR_ORACLE=ON to run the optional MPFR oracle tests in every
# matrix configuration. qa/check_oracle_matrix_cmake.sh is the convenience
# wrapper for that mode.

mkdir -p "$LOG_DIR"

exec "$CTEST_CMD" \
  -S "$SRC_DIR/qa/cmake_matrix.cmake" \
  -DSRC_DIR="$SRC_DIR" \
  -DBUILD_ROOT="$BUILD_ROOT" \
  -DLOG_DIR="$LOG_DIR"
