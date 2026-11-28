# UrQMD source discovery in CMake

When the optional UrQMD interface is enabled, CMake first checks `ISS_URQMD_SOURCE_DIR`. The default is `urqmd_code/urqmd/src` in the iSS checkout. If that directory lacks `urqmd.f`, CMake shallow-clones the Bitbucket repository into `_deps/urqmd_afterburner` under the build directory and uses its `urqmd/src` directory.

Later configurations in the same build directory reuse the clone. For an offline build or a fixed source revision, point `ISS_URQMD_SOURCE_DIR` at an existing UrQMD 3.4 checkout. A fresh clone follows the repository's default branch at configure time.

## Verification

```sh
cmake -S . -B build/urqmd-fetch-check -DISS_ENABLE_URQMD=ON \
  -DISS_URQMD_SOURCE_DIR="$PWD/build/urqmd-fetch-check/missing" \
  -DCMAKE_BUILD_TYPE=Release
cmake -S . -B build/urqmd-fetch-check -DCMAKE_DISABLE_FIND_PACKAGE_Git=ON
cmake --build build/urqmd-fetch-check -j
ctest --test-dir build/urqmd-fetch-check --output-on-failure
```

CMake cloned the Bitbucket repository at `09eeac28b5861d68d166c1f89b2e97e0a0ebfe8f` and applied the UrQMD bridge patch. The build and both UrQMD transport tests passed using that checkout. Reconfiguration succeeded with Git package discovery disabled. Separate configurations used the local checkout without fetching and skipped UrQMD discovery with `ISS_ENABLE_URQMD=OFF`. A simulated clone failure reported the error and the local-source option for offline builds.
