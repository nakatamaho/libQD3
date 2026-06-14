# QD3 consumer smoke tests

These scripts verify an installed CMake build from outside the source tree.

```sh
bash qa/consumer_smoke/run_cmake.sh /path/to/prefix
PKG_CONFIG_PATH=/path/to/prefix/lib/pkgconfig bash qa/consumer_smoke/run_pkg_config.sh /path/to/prefix
```
