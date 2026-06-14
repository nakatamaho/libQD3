#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 /path/to/install-prefix" >&2
  exit 2
fi

prefix="$1"
script_dir="$(cd "$(dirname "$0")" && pwd)"
build_dir="${BUILD_DIR:-$script_dir/build-cmake-smoke}"
rm -rf "$build_dir"
cmake -S "$script_dir" -B "$build_dir" -DCMAKE_PREFIX_PATH="$prefix"
cmake --build "$build_dir"
"$build_dir/app"
