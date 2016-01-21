\name{HBH}
\alias{HBH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Theoretical transmissions/receptions for a L-limited Hop by Hop model
}
\description{
This function computes the expected value of the number of transmissions/receptions for Hop by hop model with L-limited retransmissions per packet.
}
\usage{
HBH(p1, p2, L, N)
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
Maximum number of retransmissions.
}
  \item{N}{
Number of Hops
}
}

\details{
When there is no limitation, \code{L} value must be set as \code{L=Inf}.
}
\value{
The ouput is a matrix containing the following values for each hop and total:

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
}
\author{
Christian E. Galarza and Jonathan M. Olate
}
\seealso{
\code{\link{MCHBH}}
}
\examples{
#An N=5 Hop by hop system with limited L=7 retransmission per hop
HBH(p1=0.65,p2=0.4,L=7,N=5)

#An ilimited N=5 Hop by hop system
HBH(p1=0.65,p2=0.4,L=Inf,N=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Hop by Hop}
\keyword{transmissions}
\keyword{receptions}
\keyword{network}
