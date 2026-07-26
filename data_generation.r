
# package loading

library(Matrix) # sparse matrix
library(RSpectra) # eigen decomposition

# functions for data generation

data_generation = function(n, p, graph, dis = "gaussian", hete_ratio = 0.5, prob = 0.01, cluster_size = 50,
                           band_size = 2, con_num = 50, scale_diag = c(1, 2))
{
    # n: vector of integers, length = M, sample sizes
    # p: integer, dimension
    # graph: "diagonal", "ER", "cluster", "band"
    # dis: "gaussian" or "uniform"
    # hete_ratio: scalar in [0,1], ratio of heterogeneity
    # prob: scalar in (0,1), sparsity probability in random graph, supply if graph = "ER" or "cluster"
    # cluster_size: integer, block size of clustered graph, supply if graph = "cluster"
    # band_size: integer, width of banded graph, supply if graph = "band"
    # con_num: scalar, condition number
    # scale_diag: vector, length = 2, scale of diagonal entries of precision matrix

    M = length(n) # number of datasets

    # generate adjacency matrix (asymmetric) for "ER", "cluster", "band"

    if (graph == "ER")
    {
        A  = matrix(rbinom(p^2, 1, prob), ncol = p, nrow = p)
        diag(A) = 0
    }
    if (graph == "cluster")
    {
        num_cluster =  ceiling(p / cluster_size)
        label = rep(1:num_cluster, each = cluster_size)[1:p]

        A = matrix(0, ncol = p, nrow = p)
        for (i in 1:num_cluster)
        {
            a = which(label == i)
            A[a,a] = matrix(rbinom(length(a)^2, 1, prob), ncol = length(a), nrow = length(a))
        }
        diag(A) = 0
    }
    if (graph == "band")
    {
        A = matrix(0, ncol = p, nrow = p)
        for (i in 1:band_size)
        {
            diag(A[1:(p - i), (i + 1):p]) = 1 # upper triangle
        }
    }

    # generate partial correlation matrices

    if (graph != "diagonal")
    {
        U = matrix(runif(p^2, 1, 2), ncol = p, nrow = p)
        U = U * (matrix(rbinom(p^2, 1, 0.5), ncol = p, nrow = p) * 2 - 1)

        heter = matrix(rbinom(p^2, 1, hete_ratio), ncol = p, nrow = p) # 1 implies heterogeneity

        sigma_min = rep(0, M) # minimal eigenvalues
        sigma_max = rep(0, M) # maximal eigenvalues

        R = vector("list", length = M)
        for (m in 1:M)
        {
            U_m = matrix(runif(p^2, 1, 2), ncol = p, nrow = p)
            U_m = U_m * (matrix(rbinom(p^2, 1, 0.5), ncol = p, nrow = p) * 2 - 1)

            R_m = A * ((1 - heter) * U + heter * U_m)
            R_m = Matrix::forceSymmetric(R_m, uplo = "U")
            R_m = Matrix::Matrix(R_m, sparse = TRUE)

            sigma_max[m] = RSpectra::eigs_sym(as.matrix(R_m), k = 1, which = "LA", opts = list(retvec = FALSE))$values
            sigma_min[m] = RSpectra::eigs_sym(as.matrix(R_m), k = 1, which = "SA", opts = list(retvec = FALSE))$values

            R[[m]] = R_m
        }
        rm(A, U, heter, U_m, R_m)

        delta = (sigma_max - sigma_min) / (con_num - 1) - sigma_min
        delta_max = max(delta)
        R = lapply(R, function(R_m){
            diag(R_m) = delta_max
            R_m = R_m / delta_max
        })
    }
    if (graph == "diagonal")
    {
        R = lapply(1:M, function(m) Matrix::Matrix(diag(p), sparse = TRUE))
    }

    # generate precision matrices and datasets

    Omega_diag = runif(p, min = min(scale_diag), max = max(scale_diag))
    X = vector("list", length = M)
    Omega = vector("list", length = M)

    for (m in 1:M)
    {
        if (dis == "gaussian")
        {
            X[[m]] = matrix(rnorm(n[m] * p), nrow = n[m], ncol = p)
        }
        if (dis == "uniform")
        {
            X[[m]] = matrix(runif(n[m] * p, min = -2, max = 2), nrow = n[m], ncol = p)
        }

        Omega[[m]] = Matrix::Matrix(diag(sqrt(Omega_diag)), sparse = TRUE) %*% R[[m]] %*%
                     Matrix::Matrix(diag(sqrt(Omega_diag)), sparse = TRUE)
        temp = eigen(Omega[[m]], symmetric = TRUE)

        X[[m]] = t(temp$vectors %*% diag(1 / sqrt(temp$values)) %*% t(X[[m]]))
    }

    return(list(X = X, Omega = Omega))
}

# evaluation function for estimated precision matrices
loss_mat = function(Omega_est, Omega, sample_size, r = 1)
{
    # Omega_est: list of matrices, length = M
    # Omega: list of matrices, length = M
    # sample_size: vector, length = M
    # r: integer, >= 1

    N = sum(sample_size)

    loss_mat = Map("-", Omega_est, Omega) # list of matrices, length = M

    loss_mat_1r = Reduce("+", Map(function(x, w) w * abs(x), loss_mat, sample_size / N)) # matrix, size = p * p
    loss_mat_1r = loss_mat_1r^r

    loss_mat_2r = Reduce("+", Map(function(x, w) w * x^2, loss_mat, sample_size / N)) # matrix, size = p * p
    loss_mat_2r = loss_mat_2r^(r/2)

    loss_mat_1r = as.matrix(loss_mat_1r)
    loss_mat_2r = as.matrix(loss_mat_2r)
    p = ncol(loss_mat_1r)

    result = rep(0, 6)
    names(result) = c("loss_mat_1r_O", "loss_mat_1r_2", "loss_mat_1r_F",
                      "loss_mat_2r_O", "loss_mat_2r_2", "loss_mat_2r_F")
    result["loss_mat_1r_O"] = base::norm(loss_mat_1r, "O")
    result["loss_mat_1r_2"] = base::norm(loss_mat_1r, "2")
    result["loss_mat_1r_F"] = base::norm(loss_mat_1r, "F")^2 / p
    result["loss_mat_2r_O"] = base::norm(loss_mat_2r, "O")
    result["loss_mat_2r_2"] = base::norm(loss_mat_2r, "2")
    result["loss_mat_2r_F"] = base::norm(loss_mat_2r, "F")^2 / p

    return(result)
}
