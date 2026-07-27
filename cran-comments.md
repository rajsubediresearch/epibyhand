## Test environments

* Local: Windows 11, R 4.6.1
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1), windows-latest
  (R release), macos-latest (R release)
* win-builder: R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

## Notes

This is a new release.

The spell check flags author surnames (Haenszel, Breslow) and an
organisation acronym (IARC) in the DESCRIPTION references. All are correct
as written.

The package imports only `stats`. Scope is deliberately limited to
epidemiological measures that can be computed by hand, so that the worked
derivations the package prints correspond to a calculation a reader can
reproduce on paper.
