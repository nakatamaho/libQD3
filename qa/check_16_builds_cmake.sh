#!/usr/bin/env bash
set -u

SRC_DIR="${SRC_DIR:-$(pwd)}"
BUILD_ROOT="${BUILD_ROOT:-$SRC_DIR/_build_matrix_cmake}"
LOG_DIR="$BUILD_ROOT/logs"
CMAKE_CMD="${CMAKE:-cmake}"
CTEST_CMD="${CTEST:-ctest}"
GENERATOR="${GENERATOR:-}"
KEEP_BUILD="${KEEP_BUILD:-0}"

if [ -n "${JOBS:-}" ]; then
  jobs="$JOBS"
elif command -v nproc >/dev/null 2>&1; then
  jobs="$(nproc)"
elif command -v getconf >/dev/null 2>&1; then
  jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
else
  jobs=4
fi

mkdir -p "$BUILD_ROOT" "$LOG_DIR"
SUMMARY_OK="$BUILD_ROOT/summary.ok"
SUMMARY_NG="$BUILD_ROOT/summary.ng"
SUMMARY_TSV="$BUILD_ROOT/summary.tsv"
: > "$SUMMARY_OK"
: > "$SUMMARY_NG"
printf "tag\tresult\tlog\n" > "$SUMMARY_TSV"

fail_count=0

generator_args=()
if [ -n "$GENERATOR" ]; then
  generator_args=(-G "$GENERATOR")
fi

run_config() {
  tag="$1"
  shift
  build_dir="$BUILD_ROOT/$tag"
  log_file="$LOG_DIR/$tag.log"

  if [ "$KEEP_BUILD" != "1" ]; then
    rm -rf "$build_dir"
  fi

  echo "=== $tag ==="
  echo "build_dir: $build_dir"
  echo "log_file : $log_file"

  (
    set -e
    "$CMAKE_CMD" -S "$SRC_DIR" -B "$build_dir" "${generator_args[@]}" "$@"
    "$CMAKE_CMD" --build "$build_dir" -j"$jobs"
    "$CTEST_CMD" --test-dir "$build_dir" --output-on-failure
  ) >"$log_file" 2>&1

  rc=$?
  if [ "$rc" -eq 0 ]; then
    echo "PASS $tag"
    echo "$tag" >> "$SUMMARY_OK"
    printf "%s\tPASS\t%s\n" "$tag" "$log_file" >> "$SUMMARY_TSV"
  else
    echo "FAIL $tag"
    echo "$tag" >> "$SUMMARY_NG"
    printf "%s\tFAIL\t%s\n" "$tag" "$log_file" >> "$SUMMARY_TSV"
    tail -n 80 "$log_file" || true
    fail_count=$((fail_count + 1))
  fi
  echo
}

run_config default

for ieee_add in OFF ON; do
  for sloppy_mul in OFF ON; do
    for sloppy_div in OFF ON; do
      for fma in disabled enabled; do
        tag="ieee_add-${ieee_add}__sloppy_mul-${sloppy_mul}__sloppy_div-${sloppy_div}__fma-${fma}"
        if [ "$fma" = "enabled" ]; then
          run_config "$tag" \
            -DQD_ENABLE_IEEE_ADD="$ieee_add" \
            -DQD_ENABLE_SLOPPY_MUL="$sloppy_mul" \
            -DQD_ENABLE_SLOPPY_DIV="$sloppy_div" \
            -DQD_FMA=auto \
            -DQD_ARCH=x86-64-v3
        else
          run_config "$tag" \
            -DQD_ENABLE_IEEE_ADD="$ieee_add" \
            -DQD_ENABLE_SLOPPY_MUL="$sloppy_mul" \
            -DQD_ENABLE_SLOPPY_DIV="$sloppy_div" \
            -DQD_FMA=no \
            -DQD_ARCH=generic
        fi
      done
    done
  done
done

echo "========================================"
echo "Finished."
echo "PASS list : $SUMMARY_OK"
echo "FAIL list : $SUMMARY_NG"
echo "TSV       : $SUMMARY_TSV"
echo "Logs      : $LOG_DIR"
echo "Failures  : $fail_count"

if [ "$fail_count" -ne 0 ]; then
  exit 1
fi
