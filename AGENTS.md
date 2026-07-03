# AGENTS.md

## Scope

- Core development defaults to the legacy smoke tests (`qd_test`, `huge`,
  `pslq`, etc.) with no external test-only dependency.
- Oracle verification uses optional MPFR-backed tests under `tests/oracle`.

## Required build flow

- This branch is CMake-only. Before release-oriented testing, confirm that
  legacy generated Autotools entry points are absent:
  - `test ! -e configure && test ! -e configure.ac && test ! -e Makefile.am && test ! -e autogen.sh`

- Default smoke tests:
  - `cmake -S . -B build -DBUILD_TESTING=ON`
  - `cmake --build build -j`
  - `ctest --test-dir build --output-on-failure`

- MPFR oracle tests:
  - `cmake -S . -B build-mpfr -DBUILD_TESTING=ON -DQD3_ENABLE_MPFR_TESTS=ON`
  - `cmake --build build-mpfr -j`
  - `ctest --test-dir build-mpfr --output-on-failure`

- MPC complex oracle tests:
  - `cmake -S . -B build-mpc -DBUILD_TESTING=ON -DQD3_ENABLE_MPC_TESTS=ON`
  - `cmake --build build-mpc -j`
  - `ctest --test-dir build-mpc --output-on-failure`

## Release Gate

- The release gate must be run from a clean tree. Check:
  - `git status --short --branch`

- The MPC-inclusive release gate is the union of these CMake QA scripts:
  - `bash qa/check_16_builds_cmake.sh`
  - `bash qa/check_oracle_matrix_cmake.sh`
  - `bash qa/check_complex_oracle_matrix_cmake.sh`

- `qa/check_16_builds_cmake.sh` covers the default configuration plus the 16
  build-option combinations. It does not, by itself, enable MPFR or MPC oracle
  tests.

- `qa/check_oracle_matrix_cmake.sh` covers the same build-option matrix with
  `QD3_ENABLE_MPFR_TESTS=ON`.

- `qa/check_complex_oracle_matrix_cmake.sh` is the MPC/MPFR-backed complex
  oracle gate. It must pass for a release gate to be considered complete.

- Sanitizer and coverage runs are additional release QA. They are not a
  substitute for the MPC-inclusive release gate above.

- Sanitizer and coverage runs:
  - `cmake -S . -B build-sanitize -DQD3_ENABLE_MPFR_TESTS=ON -DQD3_ENABLE_ASAN=ON -DQD3_ENABLE_UBSAN=ON -DBUILD_TESTING=ON`
  - `cmake --build build-sanitize -j`
  - `cmake -S . -B build-coverage -DQD3_ENABLE_MPFR_TESTS=ON -DQD3_ENABLE_COVERAGE=ON -DBUILD_TESTING=ON`
  - `cmake --build build-coverage --target coverage -j`

## Seed behavior

- Oracle RNG seed resolution:
  1) `QD_TEST_SEED`
  2) `--seed=<N>`
  3) `0x9E3779B97F4A7C15`.

- Reproduction command is printed in TAP diagnostics for failures.

## Recommended QA commands

- `qa/check_16_builds_cmake.sh`
- `qa/check_oracle_matrix_cmake.sh`
- `qa/check_complex_oracle_matrix_cmake.sh`
