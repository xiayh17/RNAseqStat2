% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runDEG.R
\name{runDEG}
\alias{runDEG}
\title{run DEG analysis}
\usage{
runDEG(
  obj,
  dir = ".",
  prefix = "2-runDEG",
  parallel = TRUE,
  qc = TRUE,
  PointVolcanoParam = list(gene = 10, light = NULL, light_color = "#24ac56",
    light_label_color = "#24ac56", expend = c(0.12, 0.12))
)
}
\arguments{
\item{obj}{a DEGContainer}

\item{dir}{a directory to store results}

\item{prefix}{a prefix of file names in this step}

\item{parallel}{use parallel in DESeq2}

\item{qc}{qc in DESeq2}

\item{PointVolcanoParam}{for volcano plot}
}
\description{
run DEG analysis and summary results
}
\examples{
runDEG(DEGContainer)
}
