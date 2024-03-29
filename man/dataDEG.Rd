% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dataDEG}
\alias{dataDEG}
\title{get deg data of DEGContainer}
\usage{
dataDEG(obj, which, category = "H")
}
\arguments{
\item{obj}{a DEGContainer object}

\item{which}{kinds of DEG; can be "limma", "edgeR", "DESeq2" or "MSigDB"}

\item{category}{MSigDB collection abbreviation, such as H or C1.}
}
\value{
a dataframe
}
\description{
get deg data of DEGContainer
}
\examples{
dataDEG(DEGContainer)
}
