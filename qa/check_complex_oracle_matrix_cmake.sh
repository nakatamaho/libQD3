#!/usr/bin/env bash
set -eu

SRC_DIR="${SRC_DIR:-$(pwd)}"
BUILD_ROOT="${BUILD_ROOT:-$SRC_DIR/_build_complex_oracle_matrix_cmake}"
QD_TEST_SEED="${QD_TEST_SEED:-11400714819323198485}"
JOBS="${JOBS:-}"

for qd3_pc_dir in \
  "$HOME/mplapack/external/i/MPC/lib/pkgconfig" \
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

# The release gate requires the Fortran interfaces.  A configuration without
# a Fortran compiler must opt out explicitly with QD3_QA_FORTRAN=OFF.
QA_FORTRAN="$(printf '%s' "${QD3_QA_FORTRAN:-ON}" | tr '[:lower:]' '[:upper:]')"
case "$QA_FORTRAN" in
  ON) ;;
  OFF) echo "WARNING: QD3_QA_FORTRAN=OFF: Fortran interfaces are NOT covered by this run" >&2 ;;
  *) echo "QD3_QA_FORTRAN must be ON or OFF, got '$QA_FORTRAN'" >&2; exit 2 ;;
esac

rm -rf "$BUILD_ROOT"
mkdir -p "$BUILD_ROOT"

cmake -S "$SRC_DIR" -B "$BUILD_ROOT/default" \
  -DBUILD_TESTING=ON \
  -DQD3_ENABLE_MPC_TESTS=ON \
  -DQD_BUILD_FORTRAN="$QA_FORTRAN"

if [ -n "$JOBS" ]; then
  cmake --build "$BUILD_ROOT/default" -j "$JOBS"
else
  cmake --build "$BUILD_ROOT/default" -j
fi

QD_TEST_SEED="$QD_TEST_SEED" \
  ctest --test-dir "$BUILD_ROOT/default" --output-on-failure
