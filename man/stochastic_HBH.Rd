\name{stochastic_HBH}
\alias{stochastic_HBH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random Probabilities Monte Carlo transmissions/receptions simulations for a L-limited Hop by Hop model
}
\description{
This function compute the mean of the number of transmissions/receptions for Hop by Hop model with L-limited retransmissions per packet simulating via Monte Carlo.
}
\usage{
stochastic_HBH(dist1,p11,p12,dist2,p21,p22,L,N,M=10^5,printout=TRUE,plotspdf=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{dist1}{
For the data success probability: probability density function. Options are "uniform" and "beta".
}
  \item{p11}{
For the data success probability: lower limit of the uniform distribution (\code{dist1 == "uniform"}) or shape1 (alpha) paremHBHr of a Beta distribution (\code{dist1 == "beta"}).
}
 \item{p12}{
For the data success probability: upper limit of the uniform distribution (\code{dist1 == "uniform"}) or shape2 (beta) paremHBHr of a Beta distribution (\code{dist1 == "beta"}).
}
 \item{dist2}{
For the ACK success probability: probability density function. Options are "uniform" and "beta".
}
  \item{p21}{
For the ACK success probability: lower limit of the uniform distribution (\code{dist1 == "uniform"}) or shape1 (alpha) paremHBHr of a Beta distribution (\code{dist1 == "beta"}).
}
 \item{p22}{
For the ACK success probability: upper limit of the uniform distribution (\code{dist1 == "uniform"}) or shape2 (beta) paremHBHr of a Beta distribution (\code{dist1 == "beta"}).
}
  \item{L}{
Maximum number of retransmissions
}
  \item{N}{
Number of Hops
}
  \item{M}{
Number of Monte Carlo Simulations
}
  \item{printout}{
If \code{TRUE} (by default), the function prints some outputs and plots
}
  \item{plotspdf}{
If \code{TRUE} (by default), the function exports all plots in pdf in the working directory
}
}
\value{
The ouput is a matrix containing two elements:

\item{data}{a dataframe containing all Monte Carlo replications}
\item{stats}{descriptive statistics}
for
\item{1}{p1}
\item{2}{p2}
\item{1}{Success Probability}
\item{2}{Expected Data Transmissions}
\item{3}{Expected ACK Transmissions}
\item{4}{Expected Total Transmissions}
\item{5}{Expected Data Receptions}
\item{6}{Expected ACK Receptions}
\item{7}{Expected Total Receptions}
}
\references{
Palma, J.M.O.; Carvalho, L.D.P.; Goncalves, A.P.C.; Galarza, C.E.; De Oliveira, A.M., "Application of Control Theory Markov Systems to Minimize the Number of Transmissions in a Multi-hop Network," in Computer Aided System Engineering (APCASE), 2015 Asia-Pacific Conference on , vol., no., pp.296-301, 14-16 July 2015 <DOI:10.1109/APCASE.2015.59>

Olate, Jonathan Matias Palma, et al. "Networked control systems application: Minimization of the global number of interactions, transmissions and receptions in multi-hop network using discrete-time markovian jump linear systems." IEEE Latin America Transactions 14.6 (2016): 2675-2680.
}
\author{
Christian E. Galarza and Jonathan M. Olate
}
\seealso{
\code{\link{HBH}},\code{\link{MCHBH}}
}
\examples{
#Monte Carlo simulations for an N=5 Hop by Hop system
#with limited L=7 retransmission per hop

#We now consider p1 ~ Uniform(0.2,0.6)
dist1 = "uniform"
p11 = 0.2
p12 = 0.6

#and p2 ~ Beta(3,1)
dist2 = "beta"
p21 = 3
p22 = 1

#no outputs and plots
out = stochastic_HBH(dist1,p11,p12,dist2,p21,p22,L=7,N=5,M=5*10^3,printout=FALSE,plotspdf=FALSE)
out$data  #simulations
out$stats #resume

#uncommnet next line for outputs plots and pdf file
#out = stochastic_HBH(dist1,p11,p12,dist2,p21,p22,L=7,N=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Hop by hop}
\keyword{transmissions}
\keyword{receptions}
\keyword{network}
