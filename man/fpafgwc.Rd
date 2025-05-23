% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fpafgwc.R
\name{fpafgwc}
\alias{fpafgwc}
\title{Fuzzy Geographicaly Weighted Clustering with Flower Pollination Algorithm}
\usage{
fpafgwc(
  data,
  pop = NA,
  distmat = NA,
  ncluster = 2,
  m = 2,
  distance = "euclidean",
  order = 2,
  alpha = 0.7,
  a = 1,
  b = 1,
  error = 1e-05,
  max.iter = 100,
  randomN = 0,
  vi.dist = "uniform",
  nflow = 10,
  p = 0.8,
  gamma = 1,
  lambda = 1.5,
  delta = 0,
  ei.distr = "normal",
  flow.same = 10,
  r = 4,
  m.chaotic = 0.7,
  skew = 0,
  sca = 1
)
}
\arguments{
\item{data}{an object of data with d>1. Can be \code{matrix} or \code{data.frame}. If your data is univariate, bind it with \code{1} to get a 2 columns.}

\item{pop}{an n*1 vector contains population.}

\item{distmat}{an n*n distance matrix between regions.}

\item{ncluster}{an integer. The number of clusters.}

\item{m}{degree of fuzziness or fuzzifier. Default is 2.}

\item{distance}{the distance metric between data and centroid, the default is euclidean, see \code{\link[rdist]{cdist}} for details.}

\item{order, }{minkowski order. default is 2.}

\item{alpha}{the old membership effect with [0,1], if \code{alpha} equals 1, it will be same as fuzzy C-Means, if 0, it equals to neighborhood effect.}

\item{a}{spatial magnitude of distance. Default is 1.}

\item{b}{spatial magnitude of population. Default is 1.}

\item{error}{error tolerance. Default is 1e-5.}

\item{max.iter}{maximum iteration. Default is 500.}

\item{randomN}{random seed for initialisation (if uij or vi is NA). Default is 0.}

\item{vi.dist}{a string of centroid population distribution between \code{"uniform"} (default) and \code{"normal"}. Can be defined as \code{vi.dist=} in \code{opt_param}.}

\item{nflow}{number of flowers population. Can be defined as \code{npar=} in \code{opt_param}. Default is 10.}

\item{p}{switch probability between global and local pollination, Can be defined as \code{p} in \code{opt_param}. default is 0.8.}

\item{gamma}{Step size scaling factor. Can be defined as \code{gamma} in \code{opt_param}. Default is 1.}

\item{lambda}{Levy flights index parameter between [0,2]. Can be defined as \code{lambda} in \code{opt_param}. Default is 1.5.}

\item{delta}{Levi flights shift. Can be defined as \code{delta} in \code{opt_param}. Default is 0.}

\item{ei.distr}{distribution of random walk parameter. Can be defined as \code{ei.distr} in \code{opt_param}.}

\item{flow.same}{number of consecutive unchange to stop the iteration. Can be defined as \code{same=} in \code{opt_param}.}

\item{r}{weight in logistic chaotic between [0,4]. Can be used when \code{ei.distr='logchaotic'}. Can be defined as \code{chaos} in \code{opt_param}.}

\item{m.chaotic}{mapping parameter in kent chaotic between [0,1]. Can be used when \code{ei.distr='kentchaotic'}. Can be defined as \code{map} in \code{opt_param}.}

\item{skew}{Levy distribution skewness for random walk. Can be used when \code{ei.distr='levy'}. Can be defined as \code{skew} in \code{opt_param}.}

\item{sca}{Levy distribution scale for random walk. Can be used when \code{ei.distr='levy'}. Can be defined as \code{sca} in \code{opt_param}.}
}
\value{
an object of class \code{"fgwc"}.\cr
An \code{"fgwc"} object contains as follows:
\itemize{
\item \code{converg} - the process convergence of objective function
\item \code{f_obj} - objective function value
\item \code{membership} - membership matrix
\item \code{centroid} - centroid matrix
\item \code{validation} - validation indices (there are partition coefficient (\code{PC}), classification entropy (\code{CE}), 
SC index (\code{SC}), separation index (\code{SI}), Xie and Beni's index (\code{XB}), IFV index (\code{IFV}), and Kwon index (Kwon))
\item \code{max.iter} - Maximum iteration
\item \code{cluster} - the cluster of the data
\item \code{finaldata} - The final data (with the cluster)
\item \code{call} - the syntax called previously
\item \code{time} - computational time.
}
}
\description{
Fuzzy clustering with addition of spatial configuration of membership matrix with centroid optimization using Flower Pollination Algorithm
}
\details{
Fuzzy Geographically Weighted Clustering (FGWC) was developed by \insertCite{fgwc;textual}{naspaclust} by adding 
neighborhood effects and population to configure the membership matrix in Fuzzy C-Means. Furthermore,
the Flower Pollination Algorithm was developed by \insertCite{Yang2012;textual}{naspaclust} in order to get a more optimal
solution of a certain complex function.
}
\examples{
data('census2010')
data('census2010dist')
data('census2010pop')
# First way
res1 <- fpafgwc(census2010,census2010pop,census2010dist,3,2,'euclidean',4,nflow=10)
# Second way
# initiate parameter
param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
               alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
## tune the FPA parameter
fpa_param <- c(vi.dist='normal',npar=5,same=15,p=0.7,
               gamma=1.2,lambda=1.5,ei.distr='logchaotic',chaos=3) 
##FGWC with FPA
res2 <- fgwc(census2010,census2010pop,census2010dist,'fpa',param_fgwc,fpa_param)
}
\seealso{
\code{\link{fpafgwc}} \code{\link{gsafgwc}}
}
