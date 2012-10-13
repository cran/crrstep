crr.fast <- function (ftime, fstatus, cov1, cov2, tf, cengroup, failcode = 1, 
    cencode = 0, subset, na.action = na.omit, gtol = 1e-06, maxiter = 10, 
    init, variance=FALSE) 
{
    call <- match.call()
    cov1.name <- deparse(substitute(cov1))
    cov1.vars <- cov2.vars <- NULL
    if (!missing(cov1)) {
        cov1.vars <- colnames(as.matrix(cov1))
    }
    cov2.name <- deparse(substitute(cov2))
    if (!missing(cov2)) {
        cov2.vars <- colnames(as.matrix(cov2))
    }
    d <- data.frame(ftime = ftime, fstatus = fstatus, cengroup = if (missing(cengroup)) 
        rep(1, length(fstatus))
    else cengroup)
    if (!missing(cov1)) {
        cov1 <- as.matrix(cov1)
        nc1 <- ncol(cov1)
        d <- cbind(d, cov1)
    }
    else {
        nc1 <- 0
    }
    if (!missing(cov2)) {
        cov2 <- as.matrix(cov2)
        nc2 <- ncol(cov2)
        d <- cbind(d, cov2)
    }
    else {
        nc2 <- 0
    }
    if (!missing(subset)) 
        d <- d[subset, ]
    tmp <- nrow(d)
    d <- na.action(d)
    nmis <- 0
    if (nrow(d) != tmp) {
        nmis <- tmp - nrow(d)
        cat(format(nmis), "cases omitted due to missing values\n")
    }
    d <- d[order(d$ftime), ]
    ftime <- d$ftime
    cenind <- ifelse(d$fstatus == cencode, 1, 0)
    fstatus <- ifelse(d$fstatus == failcode, 1, 2 * (1 - cenind))
    ucg <- sort(unique.default(d$cengroup))
    cengroup <- match(d$cengroup, ucg)
    ncg <- length(ucg)
    uuu <- matrix(0, nrow = ncg, ncol = length(ftime))
    for (k in 1:ncg) {
        u <- do.call("survfit", list(formula = Surv(ftime, cenind) ~ 
            1, data = data.frame(ftime, cenind, cengroup), subset = cengroup == 
            k))
        u <- summary(u, times = sort(ftime * (1 - .Machine$double.eps)))
        uuu[k, 1:length(u$surv)] <- u$surv
    }
    uft <- sort(unique(ftime[fstatus == 1]))
    ndf <- length(uft)
    if (nc2 == 0) {
        cov1 <- as.matrix(d[, (1:nc1) + 3])
        np <- nc1
        npt <- 0
        cov2 <- 0
        tfs <- 0
    }
    else if (nc1 == 0) {
        cov2 <- as.matrix(d[, (1:nc2) + 3 + nc1])
        npt <- np <- nc2
        cov1 <- 0
        tfs <- tf(uft)
    }
    else {
        cov1 <- as.matrix(d[, (1:nc1) + 3])
        cov2 <- as.matrix(d[, (1:nc2) + 3 + nc1])
        npt <- nc2
        np <- nc1 + nc2
        tfs <- tf(uft)
    }
    if (missing(init)) 
        b <- rep(0, np)
    else b <- init
    stepf <- 0.5
    for (ll in 0:maxiter) {
        z <- .Fortran("crrfsv", as.double(ftime), as.integer(fstatus), 
            as.integer(length(ftime)), as.double(cov1), as.integer(np - 
                npt), as.integer(np), as.double(cov2), as.integer(npt), 
            as.double(tfs), as.integer(ndf), as.double(uuu), 
            as.integer(ncg), as.integer(cengroup), as.double(b), 
            double(1), double(np), double(np * np), double(np), 
            double(np), double(np * np), PACKAGE = "cmprsk")[15:17]
        if (max(abs(z[[2]]) * pmax(abs(b), 1)) < max(abs(z[[1]]), 
            1) * gtol) {
            converge <- TRUE
            break
        }
        if (ll == maxiter) {
            converge <- FALSE
            break
        }
        h <- z[[3]]
        dim(h) <- c(np, np)
        sc <- -solve(h, z[[2]])
        bn <- b + sc
        fbn <- .Fortran("crrf", as.double(ftime), as.integer(fstatus), 
            as.integer(length(ftime)), as.double(cov1), as.integer(np - 
                npt), as.integer(np), as.double(cov2), as.integer(npt), 
            as.double(tfs), as.integer(ndf), as.double(uuu), 
            as.integer(ncg), as.integer(cengroup), as.double(bn), 
            double(1), double(np), PACKAGE = "cmprsk")[[15]]
        i <- 0
        while (is.na(fbn) || fbn > z[[1]] + (1e-04) * sum(sc * 
            z[[2]])) {
            i <- i + 1
            sc <- sc * stepf
            bn <- b + sc
            fbn <- .Fortran("crrf", as.double(ftime), as.integer(fstatus), 
                as.integer(length(ftime)), as.double(cov1), as.integer(np - 
                  npt), as.integer(np), as.double(cov2), as.integer(npt), 
                as.double(tfs), as.integer(ndf), as.double(uuu), 
                as.integer(ncg), as.integer(cengroup), as.double(bn), 
                double(1), double(np), PACKAGE = "cmprsk")[[15]]
            if (i > 20) 
                break
        }
        if (i > 20) {
            converge <- FALSE
            break
        }
        b <- c(bn)
    }
    
if (variance) {
	v <- .Fortran("crrvv", as.double(ftime), as.integer(fstatus), 
        as.integer(length(ftime)), as.double(cov1), as.integer(np - 
            npt), as.integer(np), as.double(cov2), as.integer(npt), 
        as.double(tfs), as.integer(ndf), as.double(uuu), as.integer(ncg), 
        as.integer(cengroup), as.double(b), double(np * np), 
        double(np * np), double(np * np), double(length(ftime) * 
            (np + 1)), double(np), double(np), double(2 * np), 
        double(np), double(ncg * np), integer(ncg), PACKAGE = "cmprsk")[15:16]
    dim(v[[2]]) <- dim(v[[1]]) <- c(np, np)
    h <- solve(v[[1]])
    v <- h %*% v[[2]] %*% t(h)
    r <- .Fortran("crrsr", as.double(ftime), as.integer(fstatus), 
        as.integer(length(ftime)), as.double(cov1), as.integer(np - 
            npt), as.integer(np), as.double(cov2), as.integer(npt), 
        as.double(tfs), as.integer(ndf), as.double(uuu), as.integer(ncg), 
        as.integer(cengroup), as.double(b), double(ndf * np), 
        double(np), double(np), PACKAGE = "cmprsk")[[15]]
    bj <- .Fortran("crrfit", as.double(ftime), as.integer(fstatus), 
        as.integer(length(ftime)), as.double(cov1), as.integer(np - 
            npt), as.integer(np), as.double(cov2), as.integer(npt), 
        as.double(tfs), as.integer(ndf), as.double(uuu), as.integer(ncg), 
        as.integer(cengroup), as.double(b), double(ndf), double(np), 
        PACKAGE = "cmprsk")[[15]]
} else {
v <- NA
r <- NA
bj <- NA
}
    nobs <- length(ftime)
    b0 <- rep(0, length(b))
    fb0 <- .Fortran("crrf", as.double(ftime), as.integer(fstatus), 
        as.integer(length(ftime)), as.double(cov1), as.integer(np - 
            npt), as.integer(np), as.double(cov2), as.integer(npt), 
        as.double(tfs), as.integer(ndf), as.double(uuu), as.integer(ncg), 
        as.integer(cengroup), as.double(b0), double(1), double(np), 
        PACKAGE = "cmprsk")[[15]]
    if (nc1 > 0) {
        x1 <- paste(cov1.name, 1:nc1, sep = "")
        if (is.null(cov1.vars)) 
            cov1.vars <- x1
        else cov1.vars <- ifelse(cov1.vars == "", x1, cov1.vars)
    }
    if (nc2 > 0) {
        x1 <- paste(cov2.name, 1:nc2, sep = "")
        if (is.null(cov2.vars)) 
            cov2.vars <- x1
        else cov2.vars <- ifelse(cov2.vars == "", x1, cov2.vars)
        x1 <- paste("tf", 1:nc2, sep = "")
        x2 <- colnames(tfs)
        if (!is.null(x2)) 
            x1 <- ifelse(x2 == "", x1, x2)
        cov2.vars <- paste(cov2.vars, x1, sep = "*")
    }
    names(b) <- c(cov1.vars, cov2.vars)
    z <- list(coef = b, loglik = -z[[1]], score = -z[[2]], inf = matrix(z[[3]], 
        np, np), var = v, res = t(matrix(r, nrow = np)), uftime = uft, 
        bfitj = bj, tfs = as.matrix(tfs), converged = converge, 
        call = call, n = nobs, n.missing = nmis, loglik.null = -fb0)
    class(z) <- "crr"
    z
}