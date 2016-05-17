\name{edsToMatrix}
\alias{edsToMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
edsToMatrix(table, ntargets, Ctval, SampleID, miRID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
    A table from parser_of_eds.R
}
  \item{ntargets}{
     Number of targets/miRNAs each sample has
}
  \item{Ctval}{
     The column number which contains the Ctvalues~~
}
  \item{SampleID}{
     The column which contain the Sample Id
}
  \item{miRID}{
     The column which contrain the miRNA ID /name 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Anna Olofsson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (table, ntargets, Ctval, SampleID, miRID) 
{
    samp <- nrow(table)/ntargets
    data <- matrix(0, nrow = ntargets, ncol = samp)
    rownames(data) <- as.character(table[1:ntargets, miRID])
    colnames(data) <- colnames(data, do.NULL = F)
    k = ntargets
    m = 1
    for (i in seq(from = 1, to = samp)) {
        colnames(data)[i] <- table[(k * i), SampleID]
        data[, i] <- table[m:(i * k), Ctval]
        m = (i * k + 1)
    }
    data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
