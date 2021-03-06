% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linear_method.R
\name{mni152reg_mtx}
\alias{mni152reg_mtx}
\title{Get fsaverage (MNI305) to MNI152 transformation matrix.}
\usage{
mni152reg_mtx()
}
\description{
This returns the 4x4 matrix from the FreeSurfer Coordinate Systems documentation.
}
\note{
This is the opposite of using the \cite{Wu et al.} approach. It is mainly implemented in this package to allow you to easily check the difference between the methods.
}
\keyword{internal}
