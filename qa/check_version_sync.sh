#!/usr/bin/env bash
set -euo pipefail

SRC_DIR="${SRC_DIR:-$(pwd)}"
configure_ac="$SRC_DIR/configure.ac"
cmake_lists="$SRC_DIR/CMakeLists.txt"

patch="$(sed -n 's/^define(\[QD_PATCH_VERSION\], *\([0-9][0-9]*\)).*/\1/p' "$configure_ac")"
if [ -z "$patch" ]; then
  echo "ERROR: could not parse QD_PATCH_VERSION from configure.ac" >&2
  exit 1
fi

auto_base="$(sed -n 's/^AC_INIT(\[qd3\],\[\([0-9][0-9]*\)\.\([0-9][0-9]*\)\.QD_PATCH_VERSION\].*/\1.\2/p' "$configure_ac")"
if [ -z "$auto_base" ]; then
  echo "ERROR: could not parse AC_INIT qd3 version from configure.ac" >&2
  exit 1
fi

auto_version="${auto_base}.${patch}"
cmake_version="$(sed -n 's/^project(qd3 VERSION \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\) LANGUAGES.*/\1/p' "$cmake_lists")"
if [ -z "$cmake_version" ]; then
  echo "ERROR: could not parse project(qd3 VERSION ...) from CMakeLists.txt" >&2
  exit 1
fi

if [ "$auto_version" != "$cmake_version" ]; then
  echo "ERROR: version mismatch: configure.ac=$auto_version CMakeLists.txt=$cmake_version" >&2
  exit 1
fi

case "$cmake_version" in
  *.$patch) ;;
  *)
    echo "ERROR: patch version macro handling is inconsistent: patch=$patch cmake=$cmake_version" >&2
    exit 1
    ;;
esac

echo "version sync OK: $cmake_version"
