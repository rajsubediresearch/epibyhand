## Test environments

* Local: Windows 11, R 4.6.0
* GitHub Actions: ubuntu-latest (R release, R devel), windows-latest (R
  release), macos-latest (R release)
* win-builder: R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes

This is a new release.

The package imports only `stats`. Scope is deliberately limited to
epidemiological measures that can be computed by hand, so that the worked
derivations the package prints correspond to a calculation a reader can
reproduce on paper.
