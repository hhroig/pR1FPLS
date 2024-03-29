% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_seq_frob_par.R
\name{cv_seq_frob_par}
\alias{cv_seq_frob_par}
\title{Sequential crossvalidation in parallel (different penalty for each component) using Frobenius norm metric.}
\usage{
cv_seq_frob_par(
  X,
  Y,
  center = TRUE,
  argvals = NULL,
  nodes = NULL,
  nbasis = NULL,
  penalty_vec,
  ncomp = min(10, ncol(X)),
  folds = 5,
  basisobj = NULL,
  R0 = NULL,
  P = NULL,
  method = NULL,
  tol = .Machine$double.eps^0.5,
  verbose = TRUE,
  stripped = TRUE
)
}
\arguments{
\item{X}{a number-of-observations times nodes matrix.}

\item{Y}{a number-of-observations times reponses matrix.}

\item{center}{logical, indicating if data should be centered.}

\item{argvals}{a set of argument values. Not needed for FEM basis.}

\item{nodes}{a 2-column matrix with the nodes. Needed for "fpls_tps".}

\item{nbasis}{number of TPS basis to use in method "fpls_tps".}

\item{penalty_vec}{a vector of possible penalties.}

\item{ncomp}{number of components, integer.}

\item{folds}{a user defined list of folds (as generated by caret::createFolds())
or an integer indicating the number of folds.}

\item{basisobj}{a Finite Elements basis as in the fdaPDE package.}

\item{R0}{(mass) matrix of inner products between TPS basis functions.}

\item{P}{penalty matrix (optional for method "fpls_tps").}

\item{method}{only supported by: "r1fpls_fem".}

\item{tol}{convergence tolerance.}

\item{verbose}{logical, indicating if messages should be printed.}

\item{stripped}{logical.  If \code{TRUE} the calculations are stripped as
much as possible for speed. Particularly, if \code{FALSE} (default) it computes
the final models using the best combination of penalties.
Inspired by package \code{pls}.}
}
\value{
A list of crossvalidates erros (CVEs) and penalties giving the minimum
CVEs per number of components.
}
\description{
Sequential crossvalidation in parallel (different penalty for each component) using Frobenius norm metric.
}
\examples{
# 2D example:

# Generate data (30 samples, 100 nodes):
x <- seq(0, 1, length.out = 10)
y <- seq(0, 1, length.out = 10)

L <- generate_2d_data(x, y, 30, 3, 0.95)

X <- L[["X"]]
Y <- L[["Y"]]
FEM_basis <- L[["basisobj"]]

cv_fem <- cv_seq_frob_par(X = X, Y = Y, penalty_vec = c(0.01, 10, 100),
                  ncomp = 3, folds = 5, basisobj = FEM_basis,
                  method = "r1fpls_fem",
                  verbose = TRUE, stripped = FALSE)

final_fem <-cv_fem$final_model
}
