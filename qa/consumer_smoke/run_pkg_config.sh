#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 /path/to/install-prefix" >&2
  exit 2
fi

prefix="$1"
script_dir="$(cd "$(dirname "$0")" && pwd)"
build_dir="${BUILD_DIR:-$script_dir/build-pkg-config-smoke}"
mkdir -p "$build_dir"
rm -f "$build_dir/app"

if [ -d "$prefix/lib64" ]; then
  libdir="$prefix/lib64"
else
  libdir="$prefix/lib"
fi

if [ -z "${PKG_CONFIG_PATH:-}" ]; then
  export PKG_CONFIG_PATH="$libdir/pkgconfig"
fi

case "${LD_LIBRARY_PATH:-}" in
  "") export LD_LIBRARY_PATH="$libdir" ;;
  *) export LD_LIBRARY_PATH="$libdir:$LD_LIBRARY_PATH" ;;
esac

cxx="${CXX:-c++}"
$cxx "$script_dir/main.cpp" $(pkg-config --cflags --libs qd) -o "$build_dir/app"
"$build_dir/app"
