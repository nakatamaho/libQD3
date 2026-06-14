# AGENTS.md

## Scope

- Core development defaults to the legacy smoke tests (`qd_test`, `huge`,
  `pslq`, etc.) with no external test-only dependency.
- Oracle verification uses optional MPFR-backed tests under `tests/oracle`.

## Required build flow

- Default smoke tests:
  - `./configure && make check`
  - `cmake -S . -B build && cmake --build build -j`

- MPFR oracle tests (Autotools):
  - `./configure --enable-mpfr-tests`
  - `make check`

- MPFR oracle tests (CMake):
  - `cmake -S . -B build-mpfr -DQD3_ENABLE_MPFR_TESTS=ON`
  - `cmake --build build-mpfr -j`
  - `ctest --test-dir build-mpfr --output-on-failure`

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
