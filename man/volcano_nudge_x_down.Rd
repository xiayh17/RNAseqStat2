% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualization_DEGvolcano.R
\name{volcano_nudge_x_down}
\alias{volcano_nudge_x_down}
\title{nudge gene labels in x direction in down side}
\usage{
volcano_nudge_x_down(object, label_data, which, just = 0, category = "H")
}
\arguments{
\item{object}{a grouped DEGContainer}

\item{label_data}{genePoint}

\item{which}{kinds of DEG; can be "limma", "edgeR", "DESeq2" or "MSigDB"}

\item{just}{just position}

\item{category}{for \code{which} is "MSigDB"}
}
\description{
values for \code{nudge_x} in \code{\link{geom_volcano_text}}
}
