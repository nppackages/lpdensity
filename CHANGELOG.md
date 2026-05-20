# Changelog

Notable project changes are listed from newest to oldest.

## 2026-05-20 - Package Modernization

- Modernized release and repository infrastructure following the BINSREG setup, including GitHub Actions CI, Dependabot configuration, issue templates, pull request template, security policy, and Python PyPI publishing through GitHub trusted publishing.
- Renamed the active local and remote default branch from `master` to `main`, updated Stata installation URLs, and matched BINSREG GitHub repository settings for wiki, projects, issues, environments, security, and branch protection.
- Reorganized the top-level and Python READMEs, refreshed package metadata, and standardized the GPL-3 license notice.
- Updated maintainer and contributor metadata across Python, R, and Stata package files, keeping Matias D. Cattaneo as maintainer.
- Bumped package versions for the next major release: Python `3.0.0`, R `3.0`, and Stata `3.0`.
- Regenerated the Stata help PDFs for `lpdensity` and `lpbwdensity`.
- Added Python, R, and Stata numerical baseline fixtures for a deterministic local-polynomial density example.
- Added Python compatibility fixes for modern NumPy/pandas behavior in fixed-bandwidth density estimation and bandwidth-selection summaries while preserving numerical agreement with R/Stata baselines.
- Improved bandwidth-selection internals in Python, R, and Stata, including cached/reused moment calculations and a Stata ROT/IROT weight-ordering fix while preserving existing numerical baselines.
- Expanded Python, R, and Stata bandwidth-selection baselines to cover weighted DPI and ROT selectors, and reduced repeated local subsample work in Stata DPI bandwidth calculations.
