% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainfunction.R
\name{fgwc}
\alias{fgwc}
\title{Fuzzy Geographicaly Weighted Clustering}
\usage{
fgwc(data, pop, distmat, algorithm = "classic", fgwc_param, opt_param)
}
\arguments{
\item{data}{an object of data with d>1. Can be \code{matrix} or \code{data.frame}. If your data is univariate, bind it with \code{1} to get a 2 columns.}

\item{pop}{an n*1 vector contains population.}

\item{distmat}{an n*n distance matrix between regions.}

\item{algorithm}{algorithm used for FGWC}

\item{fgwc_param}{a vector that consists of FGWC parameter (see \code{\link{fgwcuv}} for parameter details)}

\item{opt_param}{a vector that consists of optimization algorithm parameter (see \code{\link{fgwcuv}} for parameter details)}
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
SC index (\code{SC}), separation index (\code{SI}), Xie and Beni's index (\code{XB}), IFV index (\code{IFV}), and Kwon index (\code{Kwon}))
\item \code{max.iter} - Maximum iteration
\item \code{cluster} - the cluster of the data
\item \code{finaldata} - The final data (with the cluster)
\item \code{call} - the syntax called previously
\item \code{time} - computational time.
}
}
\description{
Fuzzy clustering with addition of spatial configuration of membership matrix
}
\details{
Fuzzy Geographically Weighted Clustering (FGWC) was developed by \insertCite{fgwc;textual}{naspaclust} by adding 
neighborhood effects and population to configure the membership matrix in Fuzzy C-Means. There are seven optimisation algorithms that currently
provided in this package, mainly from the \insertCite{yang2014;textual}{naspaclust}. The optimization algorithm uses the centroid as the parameter to be optimized. Here are the
algorithm that can be used:
\itemize{
\item \code{"classic"} - The classical algorithm of FGWC based on \insertCite{fgwc;textual}{naspaclust} for centroid optimisation 
and \insertCite{Runkler2006;textual}{naspaclust} for membership optimization.
\item \code{"abc"} - Optimization using artificial bee colony algorithm based on \insertCite{Karaboga2007;textual}{naspaclust} 
\insertCite{@see also @fgwcabc1 and @fgwcabc2 for FGWC implementation}{naspaclust}.
\item \code{"fpa"} - Optimization using flower pollination algorithm based on \insertCite{Yang2012}{naspaclust}.
\item \code{"gsa"} - Optimization using gravitational search algorithm based on \insertCite{rashedi2009;textual}{naspaclust} and 
\insertCite{Li2017gsa;textual}{naspaclust} \insertCite{@see also @fgwcgsa for FGWC implementation}{naspaclust}.
\item \code{"hho"} - Optimization using harris-hawk optimization with \code{"heidari"} \insertCite{Heidari2019}{naspaclust} (default).
and \code{"bairathi"} \insertCite{Bairathi2018}{naspaclust}.
\item \code{"ifa"} - Optimization using intelligent firefly algorithm based on \insertCite{Yang2009;textual}{naspaclust}, 
as well as the intelligent improvement by \insertCite{intfa;textual}{naspaclust} \insertCite{@see also @Nasution2020 for FGWC implementation}{naspaclust}.
\item \code{"pso"} - Optimization using particle swarm optimization based on \insertCite{Runkler2006;textual}{naspaclust} and 
\insertCite{Bansal2011;textual}{naspaclust} for inertia option \insertCite{@see also @fgwcpso; @putra2017 for FGWC implementation}{naspaclust}.
\item \code{"tlbo"} - Optimization using teaching - learning based optimization based on \insertCite{Rao2012;textual}{naspaclust} and 
elitism improvement by \insertCite{Rao2012b;textual}{naspaclust}.
}
Furthermore, there are 10 distance that can be used to calculate the membership (see \code{\link[rdist]{cdist}} for details).
the default parameter of FGWC (in case you do not want to tune anything) is \cr \code{
c(kind='u',ncluster=2,m=2,distance='euclidean',order=2,alpha=0.7,a=1,b=1,}\cr
\code{max.iter=500,error=1e-5,randomN=1)}.\cr
There is also a universal parameter to the optimization algorithm as well as the details. The default parameter
for the optimization algorithm is \cr
\code{c(vi.dist='uniform',npar=10,par.no=2,par.dist='euclidean',par.order=2,pso=TRUE,}\cr
\code{same=10,type='sim.annealing',ei.distr='normal',vmax=0.7,wmax=0.9,wmin=0.4,}\cr
\code{chaos=4,x0='F',map=0.7,ind=1,skew=0,sca=1)} \cr
If you do not define a certain parameter, the parameter will be set to its default value
}
\examples{
data('census2010')
data('census2010dist')
data('census2010pop')
# initiate parameter
param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
               alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
## FGWC with classical algorithm
res1 <- fgwc(census2010,census2010pop,census2010dist,'classic',param_fgwc,1)
## tune the ABC parameter
abc_param <- c(vi.dist='normal',npar=5,pso=FALSE,same=15,n.onlooker=5,limit=5) 
## FGWC with ABC optimization algorithm
res2 <- fgwc(census2010,census2010pop,census2010dist,'abc',param_fgwc,abc_param) 
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{fgwcuv}}, \code{\link{abcfgwc}}, \code{\link{fpafgwc}},
\code{\link{gsafgwc}}, \code{\link{hhofgwc}}, \code{\link{ifafgwc}}, \code{\link{psofgwc}}, \code{\link{tlbofgwc}}
}
