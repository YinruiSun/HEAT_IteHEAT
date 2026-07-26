# package loading

library(Matrix) # sparse matrix
library(JGL) # JGL method
library(HGSL) # HGSL method

source("functions.r")

# ---------------------- joint graphical Lasso -------------------------------
# Danaher, P., Wang, P., & Witten, D. M. (2014). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society Series B: Statistical Methodology, 76(2), 373-397.

jgl_fun = function(X, nlambda = 50)
{
    # X: list, length = M, M sample matrices

    M = length(X)
    n = sapply(X, function(x) nrow(x)) # vector, length = M
    p = ncol(X[[1]])

    lambda2_seq = 1:nlambda / nlambda / 4

    IC = Inf
    for (i in 1:nlambda)
    {
        lambda2 = lambda2_seq[i]
        Omega_est = JGL::JGL(X, penalty = "group", 0, lambda2, weights = "sample.size", return.whole.theta = TRUE)$theta # list, length = M
        Omega_est = Map(function(x) Matrix::Matrix(x, sparse = TRUE), Omega_est)

        IC_new = 0
        for (m in 1:M)
        {
            chol_decomp = cholesky_check(Omega_est[[m]])
            if (is.null(chol_decomp) == TRUE) # check PD
            {
                IC_new = Inf
                break
            }else {
                IC_new = IC_new - 2 * n[m] * sum(log(diag(chol_decomp))) + sum(Omega_est[[m]] * cov(X[[m]])) * n[m] + log(n[m]) * (sum(Omega_est[[m]] != 0) - p) / 2
            }
        }
        if (IC_new < IC)
        {
            IC = IC_new
            lambda2_final = lambda2
            Omega_est_final = Omega_est
        }
    }

    return(list(Omega_est = Omega_est_final, lambda2 = lambda2_final))
}

# -------- heterogeneous group square-root Lasso (HGSL) ---------
# Ren, Z., Kang, Y., Fan, Y., & Lv, J. (2019). Tuning-free heterogeneous inference in massive networks. Journal of the American Statistical Association, 114(528), 1908-1925.


hgsl_fun = function(X, T_sim = 10000)
{
    # X: list of matrices, length = M, M sample matrices
    # T_sim: integer, simulation number

    M = length(X) # number of datasets
    n = sapply(X, function(x) nrow(x)) # vector, length = M
    n0 = min(n)
    N = sum(n)
    p = ncol(X[[1]])
    X = Map(function(x) scale(x, center = TRUE, scale = FALSE), X)

    grps = rep(1:(p-1), M)
    index = rep(1, 2 * M)
    index[2] = n[1]
    for (m in 2:M)
    {
        index[2*m-1] = index[2*m-2] + 1
        index[2*m] = index[2*m-2] + n[m]
    }

    lambda = rep(0, T_sim)
    for (t in 1:T_sim)
    {
        Z1 = split(rnorm(N), rep(1:M, n))
        Z2 = split(rnorm(N), rep(1:M, n))
        temp = mapply(function(z1,z2,nt) sqrt(nt) * sum(z1*z2) / sqrt( sqrt(sum(z1^2)) * sqrt(sum(z2^2)) ), Z1, Z2, n) # Z_{t,T}, t = 1,2,...,M
        lambda[t] = sqrt(sum(temp^2))
    }
    lambda = quantile(lambda, probs = 1 - 1/p, names = FALSE, type = 1)
    lambda = lambda / sqrt(n0)

    gamma_est = matrix(0, nrow = (p-1)*M, ncol = p)
    for (j in 1:p)
    {
        X_regress = as.matrix(Matrix::bdiag(Map(function(x) x[,-j], X)))
        Y_regress = Reduce(c, Map(function(x) as.vector(x[,j]), X))
        gamma_est[, j] = as.vector(HGSL::S_TISP_Path(X_regress, Y_regress, grps, M, index, lambda))
    }
    rm(X_regress, Y_regress)

    gamma_est = split(gamma_est, rep(1:M, each = p-1)) # list of vectors, length = M
    gamma_est = Map(function(x) matrix(x, nrow = p-1, ncol = p), gamma_est) # list of (p-1) by p matrix, length = M
    Omega_gamma_est = function(gamma_est)
    {
        # gamma_est: matrix, dim = (p-1) * p
        p = ncol(gamma_est)
        Omega_est = matrix(0, ncol = p, nrow = p)
        Omega_est[upper.tri(Omega_est, diag = FALSE)] = - gamma_est[upper.tri(gamma_est, diag = FALSE)]
        Omega_est[lower.tri(Omega_est, diag = FALSE)] = - gamma_est[lower.tri(gamma_est, diag = TRUE)]
        rm(gamma_est)
        diag(Omega_est) = 1
        Omega_est = Matrix::Matrix(Omega_est, sparse = TRUE)
        Omega_est
    }
    gamma_est = Map(Omega_gamma_est, gamma_est) # list of p by p matrix, length = M

    # estimating Omega with estimated diagonal entries
    Omega_diag = Map(function(x, y) 1/colMeans((x %*% y)^2), X, gamma_est) # list of p-dimensional vectors, length = M
    Omega_est = Map(function(x, y) t(t(x) * y) , gamma_est, Omega_diag)

    return(list(Omega_est = Omega_est, lambda = lambda))
}
