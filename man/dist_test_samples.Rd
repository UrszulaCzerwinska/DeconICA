% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_test_samples.R
\name{dist_test_samples}
\alias{dist_test_samples}
\title{Test impact of each Independent Component}
\usage{
dist_test_samples(A, sample.names, quant = c(0.1, 0.9), X.counts, test.type,
  thr = 0.1, isLog = NULL, return = "p.value", wide = TRUE)
}
\arguments{
\item{A}{result of \code{\link{run_fastica}} the \code{A} matrix}

\item{sample.names}{names of samples, should correspond to number of columns
of \code{A}}

\item{quant}{quantiles to use, in form of \code{c(x, y)}}

\item{X.counts}{expression data}

\item{test.type}{test of distributions to perform}

\item{thr}{threshold of maximal p.value considered 0.1 by default}

\item{isLog}{by default \code{NULL}, if \code{X} is not counts but log,
provide the base of log, for natural logarithm use \code{exp(1)}}

\item{return}{if you want to return p.values select}

\item{wide}{should the output matrix be in wide format (FALSE preferable for
plotting)}
}
\value{
returns a matrix (in long or wide) format
}
\description{
This function is applying distribution statistical test (i.e. \code{t.test},
\code{wilcox.test}) to evaluate which ICs have highest impact on differences
between samples
}
\examples{
# numerical matrix
set.seed(123)
S <- matrix(stats::rnbinom(10000, mu = 6, size = 10), 500, 80)
dat <- matrix(runif(1600,min =1, max=10 ), 80, 80, byrow = TRUE)
A <- dat / rowSums(dat)
X <- data.frame(S \%*\% A)
res_run_ica <- run_fastica(X, row.center = TRUE, n.comp = 5, overdecompose = FALSE)

#stats::t.test
dist_test_samples(A = res_run_ica$A,
sample.names = res_run_ica$samples,
X.counts = res_run_ica$log.counts,
test.type = "t.test",
isLog = 2,
return = "p.value",
thr= 0.5)

#edgeR::exactTest
dist_test_samples(A = res_run_ica$A,
sample.names = res_run_ica$samples,
X.counts = res_run_ica$log.counts,
test.type = "exactTest",
isLog = 2,
return = "p.value",
thr= 0.5)

#for plotting
res.ttest <- dist_test_samples(A = res_run_ica$A,
sample.names = res_run_ica$samples,
X.counts = res_run_ica$log.counts,
test.type = "t.test",
isLog = 2,
return = "p.value",
thr= 0.5,
wide = FALSE)

plot_dist_test(res.ttest, plot.type = "density")
plot_dist_test(res.ttest, plot.type = "line")
}
