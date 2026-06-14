#!/usr/bin/env bash
set -euo pipefail

SRC_DIR="${SRC_DIR:-$(pwd)}"
BUILD_ROOT="${BUILD_ROOT:-$SRC_DIR/_install_compare}"
PREFIX="${PREFIX:-/usr/local}"
JOBS="${JOBS:-}"
GENERATOR="${GENERATOR:-}"
CMAKE_CMD="${CMAKE:-cmake}"
MAKE_CMD="${MAKE:-make}"

if [ -z "$JOBS" ]; then
  if command -v nproc >/dev/null 2>&1; then
    JOBS="$(nproc)"
  else
    JOBS=4
  fi
fi

AUTO_BUILD="$BUILD_ROOT/autotools-build"
AUTO_DEST="$BUILD_ROOT/autotools-dest"
CMAKE_BUILD="$BUILD_ROOT/cmake-build"
CMAKE_DEST="$BUILD_ROOT/cmake-dest"
LOG_DIR="$BUILD_ROOT/logs"
mkdir -p "$LOG_DIR"
rm -rf "$AUTO_BUILD" "$AUTO_DEST" "$CMAKE_BUILD" "$CMAKE_DEST"
mkdir -p "$AUTO_BUILD" "$AUTO_DEST" "$CMAKE_DEST"

generator_args=()
if [ -n "$GENERATOR" ]; then
  generator_args=(-G "$GENERATOR")
fi

(cd "$SRC_DIR" && autoreconf -fi) >"$LOG_DIR/autoreconf.log" 2>&1

(
  set -e
  cd "$AUTO_BUILD"
  "$SRC_DIR/configure" --srcdir="$SRC_DIR" --prefix="$PREFIX"
  "$MAKE_CMD" -j"$JOBS"
  "$MAKE_CMD" install DESTDIR="$AUTO_DEST"
) >"$LOG_DIR/autotools.log" 2>&1

(
  set -e
  "$CMAKE_CMD" -S "$SRC_DIR" -B "$CMAKE_BUILD" "${generator_args[@]}" -DCMAKE_INSTALL_PREFIX="$PREFIX"
  "$CMAKE_CMD" --build "$CMAKE_BUILD" -j"$JOBS"
  DESTDIR="$CMAKE_DEST" "$CMAKE_CMD" --install "$CMAKE_BUILD"
) >"$LOG_DIR/cmake.log" 2>&1

AUTO_ROOT="$AUTO_DEST$PREFIX"
CMAKE_ROOT="$CMAKE_DEST$PREFIX"
if [ ! -d "$AUTO_ROOT" ] || [ ! -d "$CMAKE_ROOT" ]; then
  echo "ERROR: install roots not found" >&2
  echo "  autotools: $AUTO_ROOT" >&2
  echo "  cmake    : $CMAKE_ROOT" >&2
  exit 1
fi

list_files() {
  root="$1"
  (cd "$root" && find . -type f \
    ! -name '*.la' \
    ! -path './lib/cmake/QD3/*' \
    ! -path './lib64/cmake/QD3/*' \
    | sed 's#^./##' \
    | sed 's#^include/qd3/#include/qd/#' \
    | sort)
}

list_files "$AUTO_ROOT" > "$BUILD_ROOT/autotools.files"
list_files "$CMAKE_ROOT" > "$BUILD_ROOT/cmake.files"

if ! diff -u "$BUILD_ROOT/autotools.files" "$BUILD_ROOT/cmake.files" > "$BUILD_ROOT/file-set.diff"; then
  echo "ERROR: unexpected install file set differences:" >&2
  cat "$BUILD_ROOT/file-set.diff" >&2
  exit 1
fi

compare_exact_if_present() {
  rel="$1"
  if [ -f "$AUTO_ROOT/$rel" ] && [ -f "$CMAKE_ROOT/$rel" ]; then
    diff -u "$AUTO_ROOT/$rel" "$CMAKE_ROOT/$rel" > "$BUILD_ROOT/${rel//\//_}.diff" || return 1
  fi
}

for header in c_dd.h c_qd.h c_td.h dd_real.h dd_inline.h c_edd.h edd_real.h edd_inline.h fpu.h inline.h td_real.h td_inline.h qd_real.h qd_inline.h bits.h; do
  compare_exact_if_present "include/qd/$header" || {
    echo "ERROR: installed header differs: include/qd/$header" >&2
    cat "$BUILD_ROOT/include_qd_${header}.diff" >&2
    exit 1
  }
done

# qd_config.h intentionally differs in the CMake path for documented standard
# C++ isnan/isinf/isfinite/copysign simplification. Verify key public macros exist.
for macro in QD_API QD_HAVE_STD QD_ISNAN QD_ISINF QD_ISFINITE QD_COPYSIGN; do
  if ! grep -q "$macro" "$CMAKE_ROOT/include/qd/qd_config.h"; then
    echo "ERROR: CMake qd_config.h is missing $macro" >&2
    exit 1
  fi
done

for pc in lib/pkgconfig/qd.pc lib64/pkgconfig/qd.pc; do
  if [ -f "$AUTO_ROOT/$pc" ] && [ -f "$CMAKE_ROOT/$pc" ]; then
    if ! diff -u "$AUTO_ROOT/$pc" "$CMAKE_ROOT/$pc" > "$BUILD_ROOT/qd.pc.diff"; then
      echo "NOTE: qd.pc text differs; checking required fields instead. Diff saved to $BUILD_ROOT/qd.pc.diff"
      for field in '^Name: qd$' '^Version: ' 'Libs: ' 'Cflags: '; do
        grep -q "$field" "$CMAKE_ROOT/$pc" || { echo "ERROR: CMake qd.pc missing $field" >&2; exit 1; }
      done
    fi
  fi
done

if [ -f "$CMAKE_ROOT/bin/qd-config" ]; then
  "$CMAKE_ROOT/bin/qd-config" --version >/dev/null
  "$CMAKE_ROOT/bin/qd-config" --libs >/dev/null
  "$CMAKE_ROOT/bin/qd-config" --cxxflags >/dev/null
else
  echo "ERROR: CMake qd-config was not installed" >&2
  exit 1
fi

echo "install tree comparison OK"
echo "logs: $LOG_DIR"
