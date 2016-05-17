\name{norm}
\alias{norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Global mean normalization of qPCR data
}
\description{
The function subtract the mean of the vector from each element in the vector,
all values below 10 and above 38 is set to NA as these are seen as undetermined.
}
\usage{
norm(vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
  A vector of Ct values for one sample
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


## The function is currently defined as
function (vec) 
{
    vec[vec > 38] <- NA
    vec[vec < 10] <- NA
    vec - mean(na.omit(vec))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
