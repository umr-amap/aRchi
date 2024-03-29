% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/angle3d.R
\name{angle3d}
\alias{angle3d}
\title{Calculate the zenith angle from xyz coordinates}
\usage{
angle3d(o, a, b)
}
\arguments{
\item{o}{3D coordinates of the common point of the two segments}

\item{a}{3D coordinates of the other point of segment a}

\item{b}{3D coordinates of the other point of segment b}
}
\value{
The angle in degree
}
\description{
Calculate a zenith angle between two segments from 3D (i.e X,Y,Z) coordinates
}
\examples{
origin=c(0,0,0)
a=c(0,0,1)
b=c(1,0,0)

angle3d(o=origin,a=a,b=b)

}
