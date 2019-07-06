\name{MCHBH}
\alias{MCHBH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo transmissions/receptions simulations for a L-limited Hop by Hop model
}
\description{
This function compute the mean of the number of transmissions/receptions for Hop by hop model with L-limited retransmissions per packet simulating via Monte Carlo.
}
\usage{
MCHBH(p1, p2, L, N, M = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{p1}{
Data success probability
}
  \item{p2}{
ACK success probability
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
}
\value{
The ouput is a matrix containing the following values for each hop and total:

\item{1}{MC Success Probability}
\item{2}{MC Mean Data Transmissions}
\item{3}{MC Mean ACK Transmissions}
\item{4}{MC Mean Total Transmissions}
\item{5}{MC Mean Data Receptions}
\item{6}{MC Mean ACK Receptions}
\item{7}{MC Mean Total Receptions}
}
\references{
Palma, J.M.O.; Carvalho, L.D.P.; Goncalves, A.P.C.; Galarza, C.E.; De Oliveira, A.M., "Application of Control Theory Markov Systems to Minimize the Number of Transmissions in a Multi-hop Network," in Computer Aided System Engineering (APCASE), 2015 Asia-Pacific Conference on , vol., no., pp.296-301, 14-16 July 2015 <DOI:10.1109/APCASE.2015.59>

Olate, Jonathan Matias Palma, et al. "Networked control systems application: Minimization of the global number of interactions, transmissions and receptions in multi-hop network using discrete-time markovian jump linear systems." IEEE Latin America Transactions 14.6 (2016): 2675-2680.
}
\author{
Christian E. Galarza and Jonathan M. Olate
}
\seealso{
\code{\link{HBH}},\code{\link{stochastic_HBH}}
}
\examples{
#Monte Carlo simulations for an N=5 Hop by hop system
#with limited L=7 retransmission per hop

MCHBH(p1=0.65,p2=0.4,L=7,N=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Hop by Hop}
\keyword{transmissions}
\keyword{receptions}
\keyword{network}
