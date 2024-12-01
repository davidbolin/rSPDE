context("Matern rational calculations")

test_that("Precision matrix construction", {
    
    sigma <- 1
    nu <- 0.8
    range <- 0.2
    
    # create mass and stiffness matrices for a FEM discretization
    n <- 11
    x <- seq(from = 0, to = 1, length.out = n)
    
    op_cov <- matern.rational(
        loc = x, nu = nu,
        range = range, sigma = sigma, m = 2,
        parameterization = "matern"
    )
    
    # Get the precision matrix:
    tmp <- precision(op_cov, ldl = FALSE)
    cov1 <- as.vector(tmp$A%*%solve(tmp$Q, tmp$A[which.min(x),]))
    
    tmp <- precision(op_cov, ldl = TRUE)
    Q <- t(tmp$L)%*%tmp$D%*%tmp$L
    cov2 <- as.vector(tmp$A%*%solve(Q, tmp$A[which.min(x),]))
    
    reo <- sample(n,n)
    x <- x[reo]
    op_cov <- matern.rational(
        loc = x, nu = nu,
        range = range, sigma = sigma, m = 2,
        parameterization = "matern"
    )
    
    # Get the precision matrix:
    tmp <- precision(op_cov, ldl = FALSE)
    cov3 <- as.vector(tmp$A%*%solve(tmp$Q, tmp$A[which.min(x),]))
    
    tmp <- precision(op_cov, ldl = TRUE)
    Q <- t(tmp$L)%*%tmp$D%*%tmp$L
    cov4 <- as.vector(tmp$A%*%solve(Q, tmp$A[which.min(x),]))

    expect_equal(cov1, cov2, tolerance = 1e-10)
    expect_equal(cov1[reo], cov3, tolerance = 1e-10)
    expect_equal(cov1[reo], cov4, tolerance = 1e-10)

})
