# NA

## regfusionr Developer Documentation

### Unit tests and continuous integration

Note: CI status is always red because of the required package data
(about 100 MB), that is above the 5 MB limit allowed by the CRAN checks.
One needs to follow the links and manually inspect the logs to judge the
CI results atm.

[![R-CMD-check](https://github.com/dfsp-spirit/regfusionr/workflows/R-CMD-check/badge.svg)](https://github.com/dfsp-spirit/regfusionr/actions)

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/dfsp-spirit/regfusionr?branch=mainr&svg=true)](https://ci.appveyor.com/project/dfsp-spirit/regfusionr)
AppVeyor CI under Windows

## Making a new release

- Make sure all changes are logged in CHANGES file
- Bump version in DESCRIPTION
- Build package and make sure it passes CRAN tests locally. Best done
  with a recent R version, as they may have introduced even more
  annoying checks in later versions:
  `R CMD check build . && R CMD check --as-cran fsbrain_0.5.0.tar.gz`,
  or whatever version your are building
- Upload the package to
  [winbuilder](https://win-builder.r-project.org/upload.aspx) to check
  there. The service will read package metadata for your email and
  report back via mail when done.
- If everything is green both locally and on Winbuilder, submit to CRAN
  via their [package submission
  form](https://cran.r-project.org/submit.html)
- You will receive feedback from CRAN, either package was accepted or
  some version of R they test with some check still failed. Bad luck.
  You will have to modify source and do the loop again.
- Once it passes and CRAN confirms it’s on its way to the repo, tag the
  final git submit that made it into CRAN with the version,
  e.g. `git tag v0.5.0 c2hf5hjdk3` if `c2hf5hjdk3` is the commit ID.
  Check `git log --oneline` for commit IDs. When you have tagged it like
  this locally, make sure to push the tag: `git push --tags`.
- Log into github.com, and make a release there based on the tag. Copy
  relevant CHANGES section as description.

### Contributing

We accepts PRs, just create a new feature branch, the usual dance.

For larger changes, especially those requiring additional dependencies,
please get in touch before starting the actual work. It’s better to
discuss the plan first to avoid wasted efforts.
