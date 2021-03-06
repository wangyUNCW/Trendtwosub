% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{multi.freq}
\alias{multi.freq}
\title{multi.freq function}
\usage{
multi.freq(fsam)
}
\arguments{
\item{fsam}{a list in R. The order of the vectors in the list follows the order of the subsamples.}
}
\value{
count.vec it is a collection of a sequence less.count, more.count based on freq.less function.
}
\description{
This function find trend in a sample by comparing neighboring subsamples. The subsamples are stored in a list in R.
}
\details{
The first vector of data in the list will be compared with the second vector in the list by using function freq.less. Then the second vector will be compared with the 3rd vector if there is one. The statistics collected are based on computing: \deqn{\frac{1}{n_ln_{l+1}}\sum_{i=1}^{n_l}\sum_{j=1}^{n_{l+1}}1(x_{li}<x_{(l+1)j})}
}
\examples{
x1=c(1,2,4,9,0,0,NA);x2=c(1,4,9,NA);x3=c(2,5,10);
sam=list(x1,x2,x3); #
multi.freq(sam);
}
\references{
Wang, Y., Stapleton, A. E., & Chen, C. (2018). Two-sample nonparametric stochastic order inference with an application in plant physiology. Journal of Statistical Computation and Simulation, 88(14), 2668-2683.
}
