#!/bin/sh
set -eu

jobs=${JOBS:-${MAKEFLAGS:-}}
if [ -z "$jobs" ]; then
  jobs="-j2"
fi

run_case() {
  name=$1
  cppflags=$2
  builddir="buildtest-bf-$name"
  rm -rf "$builddir"
  mkdir "$builddir"
  (
    cd "$builddir"
    CPPFLAGS="$cppflags" ../configure
    make $jobs check
  )
}

run_case default ""
run_case bf "-DQD_BF"
run_case bf_add "-DQD_BF_ADD"
run_case bf_mul "-DQD_BF_MUL"
run_case bf_mul_ieee_add "-DQD_BF_MUL -DQD_IEEE_ADD"
run_case ieee_add "-DQD_IEEE_ADD"
run_case sloppy_mul "-DQD_SLOPPY_MUL"
run_case sloppy_div "-DQD_SLOPPY_DIV"
