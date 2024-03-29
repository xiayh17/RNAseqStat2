% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualization_enrichCirclize.R
\name{heatCir}
\alias{heatCir}
\title{heatmap values for annotation circlize}
\usage{
heatCir(matrix, df, look = "Count")
}
\arguments{
\item{matrix}{row names is from, column names is to}

\item{df}{a data frame contains multigroup results of enrich}

\item{look}{which column as heatmap value}
}
\value{
a list for heatmap in circlize
}
\description{
heatmap values for annotation circlize
}
\examples{
heatCir(matrix,df)
}
