#!/bin/sh
set -eu

jobs=${JOBS:-2}

run_case() {
  name=$1
  shift
  builddir="buildtest-bf-$name"
  rm -rf "$builddir"
  cmake -S . -B "$builddir" "$@"
  cmake --build "$builddir" -j "$jobs"
  ctest --test-dir "$builddir" --output-on-failure
}

run_case default
run_case bf -DQD_ENABLE_BF=ON
run_case bf_add -DQD_ENABLE_BF_ADD=ON
run_case bf_mul -DQD_ENABLE_BF_MUL=ON
run_case bf_mul_ieee_add -DQD_ENABLE_BF_MUL=ON -DQD_ENABLE_IEEE_ADD=ON
run_case ieee_add -DQD_ENABLE_IEEE_ADD=ON
run_case sloppy_mul -DQD_ENABLE_SLOPPY_MUL=ON
run_case sloppy_div -DQD_ENABLE_SLOPPY_DIV=ON
