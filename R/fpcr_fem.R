#' Functional Principal Components Regression using 2D Finite Elements basis
#'
#' @param X a number of observations times nodes matrix.
#' @param Y a number of observations times reponses matrix.
#' @param ncomp number of components, integer.
#' @param center logical, indicating if data should be centered.
#' @param basisobj a Finite Elements basis as in the \code{fdaPDE} package.
#' @param penalty_vec a vector of of penalties.
#' @param PCAvalidation validation method in \code{fdaPDE::FPCA.FEM}.
#' @param NFolds number of folds to use in \code{fdaPDE::FPCA.FEM}.
#' @param verbose logical, indicating if messages should be printed.
#' @param stripped logical.  If \code{TRUE} the calculations are stripped as
#' much as possible for speed; this is meant for use with cross-validation or
#' simulations when only the coefficients are needed.  Defaults to
#' \code{FALSE}. Inspired by package \code{pls}.
#' @param ... further arguments.  Currently not used
#'
#' @return an fpcr_fem model.
#' @export
#'
#' @examples
#' # Generate data (50 samples, 100 nodes):
#' x <- seq(0, 1, length.out = 10)
#' y <- seq(0, 1, length.out = 10)
#'
#' L <- generate_2d_data(x, y, 50, 3, 0.95)
#'
#' X <- L[["X"]]
#' Y <- L[["Y"]]
#' FEM_basis <- L[["basisobj"]]
#' pc_res <- fpcr_fem(X, Y, ncomp = 3, center = TRUE, basisobj = FEM_basis,
#'                    penalty_vec = 10^seq(-2, 4, length.out = 5),
#'                    PCAvalidation = "KFold", NFolds = 5,
#'                    stripped = FALSE )
fpcr_fem <- function(X,
                     Y,
                     ncomp = 3,
                     center = TRUE,
                     basisobj,
                     penalty_vec = c(0.0001, 0.5, 1, 10),
                     PCAvalidation = "KFold",
                     NFolds = 5,
                     verbose = TRUE,
                     stripped = FALSE,
                     ...
) {

  tictoc::tic("FPCR-FEM")

  # number of nodes and samples:
  n_samp <- nrow(X)
  n_nodes <- ncol(X)
  n_colY <- ncol(Y)

  if (center) {

    # center:
    Xc <- scale(X, scale = FALSE)
    Yc <- scale(Y, scale = FALSE)

    # get mean of Y and X to add in the end:
    Y_mean <- attr(Yc, "scaled:center")
    X_mean <- attr(Xc, "scaled:center")


  }else {

    Xc <- X
    Yc <- Y
    X_mean <- rep(0, ncol(X))
    Y_mean <- 0

  } # center data


  FPCA_solution = fdaPDE::FPCA.FEM(datamatrix = Xc,
                                   FEMbasis = basisobj,
                                   lambda = penalty_vec,
                                   validation = PCAvalidation,
                                   NFolds = NFolds,
                                   nPC = ncomp)


  # compute mass matrix:
  R0 <- mass_mat_fun(FEMbasis = basisobj)

  # coefficients for the linear regression model:
  A <- FPCA_solution[["scores"]] %*%
    Matrix::t( FPCA_solution[["loadings.FEM"]][["coeff"]] ) %*%
    R0 %*%
    FPCA_solution[["loadings.FEM"]][["coeff"]]

  # Linear model on the scores:
  data_lm <- data.frame(A = as.matrix(A), Yc = Yc)

  res_lm <- stats::lm(Yc ~ 0 + . , data = data_lm)

  # coefficients function:
  b_K <-  res_lm[["coefficients"]]
  beta_hat <- FPCA_solution[["loadings.FEM"]][["coeff"]] %*% b_K

  if (stripped) {

    ret <- list(coefficient_function = as.matrix(beta_hat),
                R0 = R0,
                Y_mean = Y_mean,
                X_mean = X_mean,
                elapsed = tictoc::toc(quiet = !verbose) )
  }else {

    # fitted values
    Y_hat <- Xc %*% R0 %*% beta_hat + Y_mean

    ret <- list(lm_model = res_lm,
                fpca_solution = FPCA_solution,
                coefficients = b_K,
                fitted.values = Y_hat,
                coefficient_function = as.matrix(beta_hat),
                R0 = R0,
                Y_mean = Y_mean,
                X_mean = X_mean,
                elapsed = tictoc::toc(quiet = !verbose)    )

  }



  class(ret) <- "fpcr_fem"

  return(ret)


} # end: FPLS function
