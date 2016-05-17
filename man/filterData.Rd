\name{filterData}
\alias{filterData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
To filter qPCR data so each miRNA is ensured to be expressed in each biogroup in more then
50\% of the samples.
}
\usage{
filterData(case, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{case}{
  A vector with Ct values for the case group, one Ct value per sample
}
  \item{control}{
  A vector with Ct values for the control group, one Ct value per sample
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
%%  ~~who you are~~
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
function (case, control) 
{
    funx <- function(row) {
        sum(is.na(row))/length(row) < 0.5
    }
    tf.mat.case <- apply(case, 1, funx)
    tf.mat.cont <- apply(control, 1, funx)
    tf.filt <- tf.mat.case & tf.mat.cont
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
