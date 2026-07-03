#!/usr/bin/env bash
set -eu

SRC_DIR="${SRC_DIR:-$(pwd)}"
BUILD_ROOT="${BUILD_ROOT:-$SRC_DIR/_build_bf_matrix_cmake}"
LOG_DIR="$BUILD_ROOT/logs"
QD_TEST_SEED="${QD_TEST_SEED:-11400714819323198485}"
JOBS="${JOBS:-}"
CMAKE_CMD="${CMAKE:-cmake}"
CTEST_CMD="${CTEST:-ctest}"
KEEP_BUILD="${KEEP_BUILD:-0}"

for qd3_pc_dir in \
  "$HOME/mplapack/external/i/MPFR/lib/pkgconfig" \
  "$HOME/mplapack/external/i/GMP/lib/pkgconfig"; do
  if [ -d "$qd3_pc_dir" ]; then
    if [ -n "${PKG_CONFIG_PATH:-}" ]; then
      PKG_CONFIG_PATH="$qd3_pc_dir:$PKG_CONFIG_PATH"
    else
      PKG_CONFIG_PATH="$qd3_pc_dir"
    fi
  fi
done
export PKG_CONFIG_PATH

if [ "$KEEP_BUILD" != "1" ]; then
  rm -rf "$BUILD_ROOT"
fi
mkdir -p "$LOG_DIR"

SUMMARY_OK="$BUILD_ROOT/summary.ok"
SUMMARY_NG="$BUILD_ROOT/summary.ng"
SUMMARY_TSV="$BUILD_ROOT/summary.tsv"
: > "$SUMMARY_OK"
: > "$SUMMARY_NG"
printf 'tag\tresult\tlog\n' > "$SUMMARY_TSV"

fail_count=0

run_case() {
  tag=$1
  shift

  build_dir="$BUILD_ROOT/$tag"
  log_file="$LOG_DIR/$tag.log"

  if [ "$KEEP_BUILD" != "1" ]; then
    rm -rf "$build_dir"
  fi
  mkdir -p "$build_dir"
  : > "$log_file"

  echo "=== $tag ==="
  echo "build_dir: $build_dir"
  echo "log_file : $log_file"

  rc=0
  {
    echo "$ $CMAKE_CMD -S $SRC_DIR -B $build_dir -DBUILD_TESTING=ON -DQD3_ENABLE_MPFR_TESTS=ON $*"
    "$CMAKE_CMD" -S "$SRC_DIR" -B "$build_dir" \
      -DBUILD_TESTING=ON \
      -DQD3_ENABLE_MPFR_TESTS=ON \
      "$@"

    if [ -n "$JOBS" ]; then
      echo "$ $CMAKE_CMD --build $build_dir -j $JOBS"
      "$CMAKE_CMD" --build "$build_dir" -j "$JOBS"
    else
      echo "$ $CMAKE_CMD --build $build_dir -j"
      "$CMAKE_CMD" --build "$build_dir" -j
    fi

    echo "$ QD_TEST_SEED=$QD_TEST_SEED $CTEST_CMD --test-dir $build_dir --output-on-failure"
    QD_TEST_SEED="$QD_TEST_SEED" \
      "$CTEST_CMD" --test-dir "$build_dir" --output-on-failure
  } >> "$log_file" 2>&1 || rc=$?

  if [ "$rc" -eq 0 ]; then
    echo "PASS $tag"
    echo "$tag" >> "$SUMMARY_OK"
    printf '%s\tPASS\t%s\n' "$tag" "$log_file" >> "$SUMMARY_TSV"
  else
    echo "FAIL $tag"
    echo "$tag" >> "$SUMMARY_NG"
    printf '%s\tFAIL\t%s\n' "$tag" "$log_file" >> "$SUMMARY_TSV"
    tail -80 "$log_file" || true
    fail_count=$((fail_count + 1))
  fi
  echo
}

run_case bf_all -DQD_ENABLE_BF=ON
run_case bf_add -DQD_ENABLE_BF_ADD=ON
run_case bf_mul -DQD_ENABLE_BF_MUL=ON
run_case bf_mul_ieee_add -DQD_ENABLE_BF_MUL=ON -DQD_ENABLE_IEEE_ADD=ON
run_case bf_all_ieee_add -DQD_ENABLE_BF=ON -DQD_ENABLE_IEEE_ADD=ON

echo "========================================"
echo "Finished."
echo "PASS list : $SUMMARY_OK"
echo "FAIL list : $SUMMARY_NG"
echo "TSV       : $SUMMARY_TSV"
echo "Logs      : $LOG_DIR"
echo "MPFR mode : ON"
echo "Seed      : $QD_TEST_SEED"
echo "Failures  : $fail_count"

if [ "$fail_count" -ne 0 ]; then
  exit 1
fi
