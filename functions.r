
# package loading

library(Matrix) # sparse matrix
library(glmnet) # lasso fit

# ------------------ functions ------------------

# t1 = Sys.time()
# output = heat(X, nite, threshold_uv = "soft", threshold_mv = "soft", splitting = FALSE)
# t2 = Sys.time()
# t2 - t1

# loss = mapply(loss_mat, output$Omega_est_integrate, list(data_gen$Omega), list(output$sample_size), 1) # matrix, dim = 6 * nite
# loss = t(loss) # matrix, dim = nite * 6
# loss_init = loss_mat(output$Omega_est_init, data_gen$Omega, output$sample_size, 1) # vector, length = 6
# loss = rbind(loss_init, loss) # matrix, dim = (nite + 1) * 6

heat = function(X, nite_max=5, ite_stop = TRUE, threshold_uv="scad", threshold_mv="scad", adapt_threshold=TRUE, nsplit=10, nlambda=50, ndelta=20, a_scad=4, a_mcp=3)
{
    # X: list, length = M, M datasets
    # nite_max: integer, maximal number of iteration
    # ite_stop: TRUE or FALSE for data-driven iteration stop when nite > 1. If FALSE, the number of iterations is nite_max; if TRUE, the number of iterations will be no more than nite_max.

    # threshold_uv: "soft", "hard", "scad" or "mcp", univariate thresholding function
    # threshold_mv: "soft", "hard", "scad" or "mcp", multivariate thresholding function
    # adapt_threshold: TRUE or FALSE, whether to estimate the variance for adaptive thresholding

    # nsplit: integer, number of data splits for cross-fitting when nite_max > 1, larger than or equal to 2
    # ndelta: integer, number of candidate shrinkage parameters in thresholding step
    # nlambda: integer, number of candidate parameters in Lasso

    # a_scad: scalar greater than 2, parameter for SCAD thresholding
    # a_mcp: scalar greater than 1, parameter for MCP thresholding

    M = length(X) # number of datasets
    sample_size = sapply(X, function(x) nrow(x)) # vector, length = M

    # individual local estimation
    local_result = Map(local_estimate, X, nlambda, adapt_threshold) # list of lists, length = M

    Omega_est_init = lapply(local_result, function(x) x$Omega_reg) # list, length = M
    Omega_debias = lapply(local_result, function(x) x$Omega_debias) # list, length = M
    IC_init = sum(sapply(local_result, function(x) x$IC)) # scalar, initial IC
    delta_lasso = sapply(local_result, function(x) x$delta_lasso) # vector, length = M
    variance_est = NULL
    if (adapt_threshold == TRUE)
    {
        variance_est = lapply(local_result, function(x) x$variance_est) # list, length = M
    }
    rm(local_result)

    # integration in central site
    center_return = center_integrate(Omega_debias, sample_size, adapt_threshold, variance_est)
    rm(variance_est)

    # shrinkage level candidates (coarse)
    Delta_uv = (1:ndelta) / ndelta * max(3, sqrt(M)) # vector, length = ndelta
    Delta_mv = (1:ndelta) / ndelta * 2 # vector, length = ndelta

    # IC in local sites
    local_IC = Map(local_shrinkage_select, X, Omega_debias, list(center_return), adapt_threshold,
                   list(Delta_uv), list(Delta_mv), threshold_uv, threshold_mv,
                   a_scad, a_mcp)

    # shrinkage level selection in central site
    delta = matrix(0, nrow = nite_max, ncol = 2) # matrix, dim = nite_max * 2
    IC = rep(0, nite_max) # vector, length = nite_max

    center_IC = center_shrinkage_select(local_IC)
    delta[1, ] = as.matrix(expand.grid(Delta_uv, Delta_mv)[center_IC$index, ])
    IC[1] = center_IC$IC

    # integrative estimation with the selected shrinkage levels
    Omega_est_integrate = vector("list", length = nite_max)
    Omega_est_integrate[[1]] = Map(local_thresholding, Omega_debias, list(center_return), adapt_threshold,
                                    delta[1, 1], delta[1, 2], threshold_uv, threshold_mv, 
                                    a_scad, a_mcp) # list, length = M
    rm(Omega_debias)

    if (nite_max > 1)
    {
        # data splitting
        fold_index = Map(function(n) sample(rep(1:nsplit, length = n)), sample_size) # list, length = M

        # cross-fitted initial estimation
        Omega_est_init_split = Map(function(fold_idx, X_mat, delta_lasso_m)
        {
            Omega_est_fold = vector("list", length = nsplit) # list, length = nsplit
            for (fold in 1:nsplit)
            {
                X_train = X_mat[which(fold_idx != fold), ]
                Omega_est_fold[[fold]] = nodewise_lasso_fixlambda(X_train, delta_lasso_m)
                # Omega_est_fold[[fold]] = nodewise_lasso(X_train, nlambda)$Omega_est
            }
            return(Omega_est_fold)
        }, fold_index, X, delta_lasso) # list, length = M

        for (t in 2:nite_max)
        {
            # iterative result in local sites
            Omega_debias = Map(local_estimate_iteration, X, Omega_est_init_split, Omega_est_integrate[[t - 1]], fold_index, nsplit) # list, length = M

            # integretion in central site
            temp = center_integrate_iteration(Omega_debias, sample_size)
            center_return$Omega_ave = temp$Omega_ave
            center_return$group_dev = temp$group_dev
            rm(temp)

            # shrinkage level candidates
            Delta_uv = delta[t - 1, 1] * 0.5 + (1:ndelta) / ndelta * delta[t - 1, 1] * 1 # vector, length = ndelta
            Delta_mv = delta[t - 1, 2] * 0.5 + (1:ndelta) / ndelta * delta[t - 1, 2] * 1 # vector, length = ndelta

            # IC in local sites
            local_IC = Map(local_shrinkage_select, X, Omega_debias, list(center_return), adapt_threshold,
                           list(Delta_uv), list(Delta_mv), threshold_uv, threshold_mv,
                           a_scad, a_mcp)  # list of vectors, length = M

            # IC selection in central site
            center_IC = center_shrinkage_select(local_IC)
            IC[t] = center_IC$IC
            delta[t, ] = as.matrix(expand.grid(Delta_uv, Delta_mv)[center_IC$index, ])

            # integrative estimation with the selected shrinkage levels
            Omega_est_integrate[[t]] = Map(local_thresholding, Omega_debias, list(center_return), adapt_threshold,
                                            delta[t, 1], delta[t, 2], threshold_uv, threshold_mv, 
                                            a_scad, a_mcp) # list, length = M

            if (ite_stop == TRUE)
            {
                # data-driven iteration stopping criterion
                if ((delta[t, 1] >= delta[t-1, 1]) & (delta[t, 2] >= delta[t-1, 2]))
                {
                    nite_max = t
                    delta = delta[1:nite_max, ]
                    IC = IC[1:nite_max]
                    Omega_est_integrate = Omega_est_integrate[1:nite_max]
                    break
                }
            }
        }
        rm(Omega_debias)
    }
    if (max(IC) == Inf) print("Warning: Iterations possibly may not converge.")

    colnames(delta) = c("delta_uv", "delta_mv")
    rownames(delta) = paste("Ite-", 1:nite_max, sep = "")
    IC = c(IC_init, IC)
    names(IC) = paste("Ite-", 0:nite_max, sep = "")

    return(list(Omega_est_integrate = Omega_est_integrate, Omega_est_init = Omega_est_init,
                delta = delta, IC = IC, nite = nite_max))
}

cholesky_check = function(A) 
{
    # cholesky decomposition & PD check
    # A: Matrix
    result = tryCatch(
    {
        Matrix::chol(A)
    },
    error = function(e)
    {
        return(NULL) # handle error
    },
    warning = function(w)
    {
        return(NULL) # handle warning
    }
    )
    return(result)
}

nodewise_lasso = function(X, nlambda)
{
    # X: matrix, dim = n * p
    # nlambda: integer, number of candidate parameters in Lasso

    p = ncol(X)
    n = nrow(X)
    X = scale(X, center = TRUE, scale = FALSE)

    # lasso fitting
    lambda_seq = 0.5 + (nlambda : 1) / nlambda * 2 # vector, length = nlambda, decreasing sequence. small values (<0.5) will not be considered

    gamma_seq = matrix(0, nrow = (p - 1) * p, ncol = nlambda)
    Omega_diag_seq = matrix(0, nrow = p, ncol = nlambda)

    # nodewise lasso procedure with given lambda sequence
    for (j in 1:p)
    {
        lambda_j = sqrt(sum(X[, j]^2) / n * log(p) / n) * lambda_seq
        fit_j = glmnet::glmnet(x = X[, -j], y = X[, j], family = "gaussian", lambda = lambda_j, intercept = FALSE)
        # fit_j$beta: matrix, dim = (p - 1) * nlambda
        gamma_seq[((j-1)*(p-1)+1):(j*(p-1)), ] = as.matrix(fit_j$beta)
        Omega_diag_seq[j, ] = n / as.vector(crossprod(X[, j] - predict(fit_j, newx = X[, -j]), X[, j])) # vector, length = nlambda
    }
    rm(lambda_j, fit_j)
    gamma_seq = Matrix::Matrix(gamma_seq, sparse = TRUE)

    # lambda selection
    IC = Inf

    for (i in 1:nlambda)
    {
        Omega = diag(p)
        gamma_temp = matrix(gamma_seq[, i], nrow = (p - 1), ncol = p)

        is_off_diag = (row(Omega) != col(Omega))
        Omega[is_off_diag] = - gamma_temp
        Omega = t(t(Omega) * Omega_diag_seq[, i])
        # symmetrization
        Omega_sym = Omega * (abs(Omega) <= abs(t(Omega))) + t(Omega) * (abs(Omega) > abs(t(Omega)))

        Omega = Matrix::Matrix(Omega, sparse = TRUE)
        Omega_sym = Matrix::Matrix(Omega_sym, sparse = TRUE)

        if (i == 1)
        {
            Omega_est = Omega
            index = 1
        }

        # IC
        chol_decomp = cholesky_check(Omega_sym)
        if (is.null(chol_decomp) == FALSE) # check PD
        {
            IC_new = 2 * sum(log(diag(chol_decomp)))
            IC_new = - n * IC_new + sum(Omega_sym * cov(X)) * n + log(n) * (sum(Omega_sym != 0) - p) / 2 # BIC
            if (IC_new <= IC)
            {
                IC = IC_new
                index = i
                Omega_est = Omega
            }
        }
    }
    rm(gamma_seq, Omega_diag_seq, Omega, Omega_sym, gamma_temp)

    return(list(Omega_est = Omega_est, delta_lasso = lambda_seq[index], IC = IC))
}

nodewise_lasso_fixlambda = function(X, lambda)
{
    # X: matrix, dim = n * p
    # lambda: scalar

    p = ncol(X)
    n = nrow(X)
    X = scale(X, center = TRUE, scale = FALSE)

    gamma = matrix(0, nrow = (p - 1), ncol = p)
    Omega_diag = rep(0, p)

    for (j in 1:p)
    {
        lambda_j = sqrt(sum(X[, j]^2) / n * log(p) / n) * lambda
        fit_j = glmnet::glmnet(x = X[, -j], y = X[, j], family = "gaussian", lambda = lambda_j, intercept = FALSE)
        gamma[, j] = as.vector(fit_j$beta)
        Omega_diag[j] = n / sum((X[, j] - predict(fit_j, newx = X[, -j])) * X[, j])
    }
    rm(fit_j)

    Omega_est = diag(p)
    is_off_diag = (row(Omega_est) != col(Omega_est))
    Omega_est[is_off_diag] = - gamma
    rm(gamma, is_off_diag)

    Omega_est = Matrix::Matrix(Omega_est, sparse = TRUE)
    Omega_est = t(t(Omega_est) * Omega_diag)

    return(Omega_est)
}

local_estimate = function(X, nlambda, adapt_threshold)
{
    # X: matrix, dim = n * p
    # nlambda: integer, number of candidate parameters in Lasso
    # adapt_threshold: TRUE or FALSE, adaptive thresholding or not. If TRUE, variance estimation will be conducted.

    n = nrow(X)
    p = ncol(X)

    # Omega estimation
    result = nodewise_lasso(X, nlambda)
    Omega_reg = result$Omega_est
    delta_lasso = result$delta_lasso
    IC = result$IC
    rm(result)

    gamma = t(t(Omega_reg) / diag(Omega_reg)) # matrix, dim = p * p, diag = 1, off-diag = -gamma
    gamma = Matrix::Matrix(gamma, sparse = TRUE)

    # Omega debiasing
    Omega_debias = Omega_reg + t(Omega_reg) - t(Omega_reg) %*% ((n - 1) / n * cov(X)) %*% Omega_reg
    Omega_debias = as.matrix(Omega_debias)
    Omega_debias = (Omega_debias + t(Omega_debias)) / 2 # mathematical equivalence

    variance_est = NULL
    if (adapt_threshold == TRUE)
    {
        variance_est = matrix(0, ncol = p, nrow = p)

        # variance estimation for adaptive thresholding
        Omega_est_temp = t(t(gamma) * diag(Omega_debias))
        Omega_epsilon_est = scale(X, center = TRUE, scale = FALSE) %*% Omega_est_temp # matrix, dim = n * p
        rm(gamma, Omega_est_temp)

        # avoid the use of apply function
        for (i in 1:n)
        {
            variance_est = variance_est + (tcrossprod(Omega_epsilon_est[i, ]) - Omega_debias)^2
        }
        rm(Omega_epsilon_est)
        variance_est = variance_est / n
    }

    return(list(Omega_reg = Omega_reg, Omega_debias = Omega_debias, delta_lasso = delta_lasso, IC = IC, variance_est = variance_est))
}

center_integrate = function(Omega_debias, sample_size, adapt_threshold, variance_est=NULL)
{
    # Omega_debias: list, length = M
    # sample_size: vector, length = M
    # adapt_threshold: TRUE or FALSE, adaptive thresholding or not
    # variance_est: list, length = M, variance estimation, supply when adapt_threshold = TRUE

    N = sum(sample_size) # total sample size
    p = ncol(Omega_debias[[1]])
    M = length(sample_size)

    # averaging
    Omega_ave = Reduce("+", Map("*", Omega_debias, sample_size / N))

    # group deviation
    group_dev = Reduce("+", Map(function(x, w) w * (x - Omega_ave)^2, Omega_debias, sample_size / N))
    group_dev = sqrt(group_dev)

    # variance aggregate
    if (adapt_threshold == TRUE)
    {
        # adaptive thresholding
        var_w_l1 = Reduce("+", Map("*", variance_est, sample_size / N)) # matrix, size = p * p
        var_l1 = Reduce("+", variance_est) # matrix, size = p * p
        var_l2 = sqrt(Reduce("+", Map(function(x) x^2, variance_est))) # matrix, size = p * p
        var_linfty = Reduce(pmax, variance_est) # matrix, size = p * p
        rm(variance_est)
    }else
    {
        # universal thresholding
        # rough surrogates to determine the scale of thresholding levels
        variance_est = lapply(Omega_debias, function(x) x^2 + tcrossprod(diag(x))) # surrogate variance estimation

        var_linfty = Reduce(pmax, variance_est) # matrix, size = p * p
        var_linfty = mean(var_linfty) # scalar
        rm(variance_est)

        var_w_l1 = NULL
        var_l1 = NULL
        var_l2 = NULL
    }

    return(list(N = N, M = M, Omega_ave = Omega_ave, group_dev = group_dev,
                var_w_l1 = var_w_l1, var_l1 = var_l1,
                var_l2 = var_l2, var_linfty = var_linfty))
}

local_estimate_iteration = function(X, Omega_est_init_split, Omega_est_integrate, fold_index, nsplit)
{
    # X: matrix, size = n * p
    # Omega_est_init_split: list, length = nsplit, cross-fitted initial estimation
    # Omega_est_integrate: matrix, size = p * p, integrative estimation in the previous iteration
    # fold_index: vector, length = n, data split index
    # nsplit: integer, number of data splits

    n = nrow(X)
    p = ncol(X)
    n_fold = sapply(1:nsplit, function(fold) sum(fold_index == fold)) # vector, length = nsplit

    Omega_debias = matrix(0, nrow=p, ncol=p)
    for (fold in 1:nsplit)
    {
        Omega_debias_temp = t(Omega_est_init_split[[fold]]) + Omega_est_integrate - t(Omega_est_init_split[[fold]]) %*% ((n_fold[fold] - 1) / n_fold[fold] * cov(X[which(fold_index == fold), ])) %*% Omega_est_integrate
        Omega_debias_temp = as.matrix(Omega_debias_temp)
        Omega_debias = Omega_debias + n_fold[fold] / n * Omega_debias_temp
    }
    rm(Omega_debias_temp)
    Omega_debias = (Omega_debias + t(Omega_debias)) / 2

    return(Omega_debias)
}

center_integrate_iteration = function(Omega_debias, sample_size)
{
    # Omega_debias: list, length = M
    # sample_size: vector, length = M

    N = sum(sample_size)

    # averaging
    Omega_ave = Reduce("+", Map("*", Omega_debias, sample_size / N))

    # group deviation
    group_dev = Reduce("+", Map(function(x, w) w * (x - Omega_ave)^2, Omega_debias, sample_size / N))
    group_dev = sqrt(group_dev)

    return(list(Omega_ave = Omega_ave, group_dev = group_dev))
}

local_thresholding = function(Omega_debias, center_return, adapt_threshold, delta_uv, delta_mv, threshold_uv, threshold_mv, a_scad, a_mcp)
{
    # Omega_debias: matrix, result from "local_estimate"
    # center_return: list, returns from "center_integrate"
    # adapt_threshold: TRUE or FALSE, adaptive thresholding or not

    # delta_uv, delta_mv: scalar, shrinkage levels for thresholding functions T_uv and T_mv
    # threshold_uv: "soft", "hard", "scad" or "mcp" for univariate thresholding
    # threshold_mv: "soft", "hard", "scad" or "mcp" for multivariate thresholding

    # a_scad: scalar greater than 2, parameter for SCAD thresholding
    # a_mcp: scalar greater than 1, parameter for MCP thresholding

    p = ncol(Omega_debias)

    Omega_ave = center_return$Omega_av
    group_dev = center_return$group_dev

    if (adapt_threshold == TRUE)
    {
        lambda_uv = delta_uv * sqrt(center_return$var_w_l1 * log(p) / center_return$N) # matrix (dim = p * p)

        lambda_mv = center_return$var_l1 + 2 * sqrt(2) * center_return$var_l2 * sqrt(log(p)) +
                    4 * center_return$var_linfty * log(p) # matrix (dim = p * p)
        lambda_mv = delta_mv * sqrt(lambda_mv / center_return$N)
    }
    if (adapt_threshold == FALSE)
    {
        lambda_uv = delta_uv * sqrt(center_return$var_linfty * log(p) / center_return$N) # scalar

        lambda_mv = center_return$var_linfty * (center_return$M + log(p)) / center_return$N # scalar
        lambda_mv = delta_mv * sqrt(lambda_mv)
    }

    # univariate thresholding
    if (threshold_uv == "soft")
    {
        Gamma_est = Omega_ave * pmax(1 - lambda_uv / abs(Omega_ave), 0)
    }
    if (threshold_uv == "hard")
    {
        Gamma_est = Omega_ave * (abs(Omega_ave) > lambda_uv)
    }
    if (threshold_uv == "scad")
    {
        Gamma_est = Omega_ave * pmax(1 - lambda_uv / abs(Omega_ave), 0) * (abs(Omega_ave) <= 2 * lambda_uv) +
                    ((a_scad - 1) * Omega_ave - a_scad * sign(Omega_ave) * lambda_uv) / (a_scad - 2) * ((2 * lambda_uv < abs(Omega_ave))&(abs(Omega_ave) <= a_scad * lambda_uv)) + 
                    Omega_ave * (abs(Omega_ave) > a_scad * lambda_uv)
    }
    if (threshold_uv == "mcp")
    {
        Gamma_est = Omega_ave * pmax(1 - lambda_uv / abs(Omega_ave), 0) * a_mcp / (a_mcp - 1) * (abs(Omega_ave) <= a_mcp * lambda_uv) + 
                    Omega_ave * (abs(Omega_ave) > a_mcp * lambda_uv)
    }
    Gamma_est = Matrix::Matrix(Gamma_est, sparse = TRUE)

    # multivariate thresholding
    if (threshold_mv == "soft")
    {
        Lambda_est = (Omega_debias - Omega_ave) * pmax(1 - lambda_mv / group_dev, 0)
    }
    if (threshold_mv == "hard")
    {
        Lambda_est = (Omega_debias - Omega_ave) * (group_dev > lambda_mv)
    }
    if (threshold_mv == "scad")
    {
        Lambda_est = (Omega_debias - Omega_ave) * pmax(1 - lambda_mv / group_dev, 0) * (group_dev <= 2 * lambda_mv) +
                     (Omega_debias - Omega_ave) * ((a_scad - 1)  - a_scad * lambda_mv / group_dev) / (a_scad - 2) * ((2 * lambda_mv < group_dev)&(group_dev <= a_scad * lambda_mv)) +
                     (Omega_debias - Omega_ave) * (group_dev > a_scad * lambda_mv)
        loc = which(group_dev == 0)
        Lambda_est[loc] = 0 # calibrate the possible cases -Inf * 0 = NaN
    }
    if (threshold_mv == "mcp")
    {
        Lambda_est = (Omega_debias - Omega_ave) * pmax(1 - lambda_mv / group_dev, 0) * a_mcp / (a_mcp - 1) * (group_dev <= a_mcp * lambda_mv) +
                     (Omega_debias - Omega_ave) * (group_dev > a_mcp * lambda_mv)
    }
    Lambda_est = Matrix::Matrix(Lambda_est, sparse = TRUE)

    Omega_est_integrate = Gamma_est + Lambda_est
    return(Omega_est_integrate)
    # return(list(Gamma_est = Gamma_est, Lambda_est = Lambda_est))
}

local_shrinkage_select = function(X, Omega_debias, center_return, adapt_threshold, Delta_uv, Delta_mv, threshold_uv, threshold_mv, a_scad, a_mcp)
{
    # X: matrix, dim = n * p

    # Omega_debias: matrix, result from "local_est_debias"
    # center_return: list, returns from "center_integrate"
    # adapt_threshold: TRUE or FALSE, adaptive thresholding or not

    # Delta_uv: sequence of candidate shrinkage parameters
    # Delta_mv: sequence of candidate shrinkage parameters

    # threshold_uv: "soft", "hard", "scad" or "mcp" for univariate thresholding
    # threshold_mv: "soft", "hard", "scad" or "mcp" for multivariate thresholding

    # a_scad: scalar greater than 2, parameter for SCAD thresholding
    # a_mcp: scalar greater than 1, parameter for MCP thresholding

    n = nrow(X)
    p = ncol(X)

    delta_grid = expand.grid(Delta_uv, Delta_mv) # matrix-like, ncol = 2
    local_IC = rep(0, nrow(delta_grid)) # vector, length = length(Delta_uv) * length(Delta_mv)

    for (i in 1:nrow(delta_grid))
    {
        delta_uv = delta_grid[i, 1]
        delta_mv = delta_grid[i, 2]

        Omega_est_integrate = local_thresholding(Omega_debias, center_return, adapt_threshold, delta_uv, delta_mv, threshold_uv, threshold_mv, a_scad, a_mcp)

        chol_decomp = cholesky_check(Omega_est_integrate)
        # IC
        if (is.null(chol_decomp) == TRUE) # check PD
        {
            local_IC[i] = Inf
        }else{
            local_IC[i] = 2 * sum(log(diag(chol_decomp)))
            local_IC[i] = - n * local_IC[i] + sum(Omega_est_integrate * cov(X)) * n
            local_IC[i] = local_IC[i] + log(n) * (sum(Omega_est_integrate != 0) - p) / 2
        }
    }

    return(local_IC)
}

center_shrinkage_select = function(local_IC)
{
    # local_IC: list of vectors, length = M

    IC = Reduce("+", local_IC)
    minimum = min(IC, na.rm = TRUE)
    if (minimum == Inf)
    {
        index = length(IC)
    }else
    {
        # index = tail(which(IC == minimum), 1) # tail: avoid the first
        index = which.min(IC)
    }
    IC = IC[index] # scalar

    return(list(index = index, IC = IC))
}

