

### Mandatory:
# Weights = NxN similarity matrix
# Method = "NJW", "Ncut", "DiffusionMap"


### Optional:
# kmeans.method="kmeans" (default)
#               "fuzzy"
#               "poly.fuzzy"
# m = fuzziness parameter
    # defaults to m=2 for fuzzy, m=.5 for poly.fuzzy

### Output:
# outputs the kmeans object; "cluster" ID's object contained somewhere 
#      within, depending on choice of kmeans.method
# => use ls(_) on the output to see all objects contained within the output


clustering <- function(Weights, method, k, t=NULL, sparse=T, 
    kmeans.method="kmeans", m=NULL) {
    
    a <- Sys.time()
    
    require(Matrix)
    
    if (method == "NJW") {
        
        D <- rowSums(Weights)   # Degrees matrix
        
        ### 
        ### Check: if row has zero similarity, problems arise
        ###
        
        if( min(D) <= 0) { 
            
            resp <- readline(prompt="One of your similarity rows has zero weight. Would you like to set 
                a 1 on the diagonal of the similarity? Type Y or N \n")
            if (resp == "Y" | resp == "y") {
                n <- which(D == 0)
                D[n] <- 1 
            }
            
            else { stop("One of your rows has zero weight.") }
        }
        
        ###
        ### Subspace Projection
        ###
        
        D <- Diagonal(n=nrow(Weights),(D^-.5))
        Z <- D %*% Weights %*% D
        
        
        if (kmeans.method=="RKM") {
            ###
            ### Run RKM on the normalized Z matrix
            ###
            
            cluster.out = RKM(Z, nclus=k, ndim=k, method = "RKM", rotation = "varimax", nstart=10)
            b <- Sys.time()
            time.elapsed <- b - a
            print(time.elapsed)
            
            return(cluster.out)
        }
        
        
        # RSpectra is efficient for dense matrices
        
        if (sparse==F) {
            require(RSpectra)
            EZ <- eigs_sym(Z, k+1, 'LM')$vector 
        }
        
        
        # irlba is efficient and accurate for sparse matrices
        
        else {
            require(irlba)
            EZ <- partial_eigen(x=Z, n = k+1, symmetric = TRUE)$vectors
        }
        
        
        # U is the L2-normalized eigenspace
        
        U <- EZ/sqrt(rowSums(EZ^2))
        
        
        ###
        ### k-Mmeans in this normalized eigenspace:
        ###
        
        # Regular k-means
        if (kmeans.method=="kmeans") {
            cluster.out <- kmeans(U, centers=k, nstart = 100)
        }
        
        # Fuzzy k-means
        else if (kmeans.method=="fuzzy") {
            require(fclust)
            if (is.null(m)) {
                m <- 2
            }
            cluster.out <- FKM(X=U,k=k, m=m, RS=10)
        }
        
        # Polynomial fuzzy k-means
        else if (kmeans.method=="poly.fuzzy") {
            require(fclust)
            if (is.null(m)) {
                m <- .5
            }
            cluster.out <- FKM.pf(X=U,k=k, b=m, RS=10)
        }
        
        else {stop("Pick a valid k-means method.")}
        
        
    }
    
    
    else if (method == "Ncut") {
        
        n <- nrow(Weights)
        dvec_inv = 1/sqrt(rowSums(Weights))
        #W_tilde = Matrix(rep(dvec_inv,n), ncol=n) * Weights * t(Matrix(rep(dvec_inv,n),ncol=n))
        W_tilde = Diagonal(n,dvec_inv) %*% Weights %*% Diagonal(n,dvec_inv)
        W_tilde = (W_tilde+t(W_tilde))/2
        
        # diag(dvec_inv) %*% Weights %*% diag(dvec_inv) ?
        # why the average part?
        
        if (sparse==F) {
            require(RSpectra)
            EZ <- eigs_sym(W_tilde, k, 'LM')$vector }
        else {
            require(irlba)
            EZ <- partial_eigen(x=W_tilde, n = k, symmetric = TRUE)$vectors
        }
        
        V <- EZ
        V = matrix(rep(dvec_inv,k-1), ncol = k-1) * V[,2:k]
        V = V / (matrix(rep(sqrt(rowSums(V^2)),k-1),ncol=k-1))
        
        
        if (kmeans.method=="kmeans") {
            cluster.out <- kmeans(V, centers=k, nstart = 100)
        }
        
        else if (kmeans.method=="fuzzy") {
            require(fclust)
            if (is.null(m)) {
                m <- 2
            }
            cluster.out <- FKM(X=V,k=k, m=m, RS=10)
        }
        
        else if (kmeans.method=="poly.fuzzy") {
            require(fclust)
            if (is.null(m)) {
                m <- .5
            }
            cluster.out <- FKM.pf(X=V,k=k, b=m, RS=10)
        }
        
        else {stop("Pick a valid k-means method.")}
        
    }
    
    
    else if (method == 'DiffusionMap'){
        
        if (is.null(t)) {
            
            stop("Specify a t value.") }
        
        require(RSpectra)
        
        n <- nrow(Weights)
        dvec_inv = 1/sqrt(rowSums(Weights))
        #W_tilde = matrix(rep(dvec_inv,n), ncol=n) * Weights * t(matrix(rep(dvec_inv,n),ncol=n))
        W_tilde = Diagonal(n,dvec_inv) %*% Weights %*% Diagonal(n,dvec_inv)
        W_tilde = (W_tilde+t(W_tilde))/2
        
        
        # diag(dvec_inv) %*% Weights %*% diag(dvec_inv) ?
        # why the average part?
        
        if (kmeans.method=="RKM") {
            ###
            ### Run RKM on the normalized W_tilde matrix
            ###
            
            cluster.out = RKM(W_tilde, nclus=k, ndim=k+1, method = "RKM", rotation = "varimax", nstart=10)
            b <- Sys.time()
            time.elapsed <- b - a
            print(time.elapsed)
            
            return(cluster.out)
        }
        
        EV <- eigs_sym(W_tilde, k+1, 'LM')
        V <- EV$vector
        lambda <- EV$value
        V_inv = 1/sqrt((rowSums(V[,2:(k+1)]^2)))
        V <- matrix(rep(V_inv,k), ncol=k) * V[,2:(k+1)]
        V = matrix(rep(dvec_inv,k), ncol = k)  * V
        V = (matrix(rep(lambda[2:(k+1)], each=n), ncol=k)^t )* V
        
        # run kmeans in eigenspace:
        
        if (kmeans.method=="kmeans") {
            cluster.out <- kmeans(V, centers=k, nstart = 100)
        }
        
        else if (kmeans.method=="fuzzy") {
            require(fclust)
            if (is.null(m)) {
                m <- 2
            }
            cluster.out <- FKM(X=V,k=k, m=m, RS=10)
        }
        
        else if (kmeans.method=="poly.fuzzy") {
            require(fclust)
            if (is.null(m)) {
                m <- .5
            }
            cluster.out <- FKM.pf(X=V,k=k, b=m, RS=10)
        }
        
        else if (kmeans.method=="RKM") {
            require(clustrd)
            cluster.out = cluspca(V, nclus=k, ndim=k, method = "RKM", rotation = "varimax", nstart=10)
        }
        
    }
    
    else {
        stop("Pick a valid clustering method.") }
        
    b <- Sys.time()
    time.elapsed <- b - a
    print(time.elapsed)
    
    return(cluster.out) 
    
}

    
    
    
    
###
### Define the RKM function
###

RKM <- function (data, nclus, ndim, alpha = NULL, method = "RKM", center = TRUE, 
    scale = TRUE, rotation = "none", nstart = 10, smartStart = NULL, 
    seed = 1234) {

    require(ggplot2)
    require(dummies)
    require(grid)
    require(corpcor)

    ssq = function(a) {
        t(as.vector(c(as.matrix(a))))%*%as.vector(c(as.matrix(a)))
    }

    if (is.null(alpha) == TRUE) {
        if (method == "RKM") {
            alpha = 0.5
        }
        else if (method == "FKM") {
            alpha = 0
        }
    }
    odata = data
    data = scale(data, center = center, scale = scale)
    # data = data.matrix(data)
    n = dim(data)[1]
    m = dim(data)[2]
    conv = 1e-06
    func = {
    }
    index = {
    }
    AA = {
    }
    FF = {
    }
    YY = {
    }
    UU = {
    }

    require(irlba)

    for (run in c(1:nstart)) {
        if (is.null(smartStart)) {
            myseed = seed + run
            set.seed(myseed)
            randVec = matrix(ceiling(runif(n) * nclus), n, 1)
        }
        else {
            randVec = smartStart
        }
        U = dummy(randVec)
        P = U %*% pseudoinverse(t(U) %*% U) %*% t(U)


        # A = eigen(t(data) %*% ((1 - alpha) * P - (1 - 2 * alpha) * 
        #         diag(n)) %*% data)$vectors
        # A = A[, 1:ndim]


        testobj <- t(data) %*% ((1 - alpha) * P - (1 - 2 * alpha) * diag(n)) %*% data

        A <- partial_eigen(x=testobj, n = ndim, symmetric = TRUE)$vectors


        G = data %*% A
        Y = pseudoinverse(t(U) %*% U) %*% t(U) %*% G
        f = alpha * ssq(data - G %*% t(A)) + (1 - alpha) * ssq(data %*% 
                A - U %*% Y)
        f = as.numeric(f)
        fold = f + 2 * conv * f
        iter = 0
        while (f < fold - conv * f) {
            fold = f
            iter = iter + 1
            outK = try(kmeans(G, centers = Y, nstart = 100), 
                silent = T)
            if (is.list(outK) == FALSE) {
                outK = EmptyKmeans(G, centers = Y)
            }
            v = as.factor(outK$cluster)
            U = diag(nlevels(v))[v, ]
            P = U %*% pseudoinverse(t(U) %*% U) %*% t(U)

            # A = eigen(t(data) %*% ((1 - alpha) * P - (1 - 2 * 
            #         alpha) * diag(n)) %*% data)$vectors
            # A = A[, c(1:ndim)]

            testobj <- t(data) %*% ((1 - alpha) * P - (1 - 2 * alpha) * diag(n)) %*% data

            A <- partial_eigen(x=testobj, n = ndim, symmetric = TRUE)$vectors


            G = data %*% A
            Y = pseudoinverse(t(U) %*% U) %*% t(U) %*% G
            f = alpha * ssq(data - G %*% t(A)) + (1 - alpha) * 
                ssq(data %*% A - U %*% Y)
        }
        func[run] = f
        FF[[run]] = G
        AA[[run]] = A
        YY[[run]] = Y
        UU[[run]] = U
        cat("Just finished iteration ", run, "\n")
    }
    mi = which.min(func)
    U = UU[[mi]]
    cluID = apply(U, 1, which.max)
    csize = round((table(cluID)/sum(table(cluID))) * 100, digits = 2)
    aa = sort(csize, decreasing = TRUE)
    require(plyr)
    cluID = mapvalues(cluID, from = as.integer(names(aa)), to = as.integer(names(table(cluID))))
    centroid = YY[[mi]]
    centroid = centroid[as.integer(names(aa)), ]
    if (rotation == "varimax") {
        require(stats)
        AA[[mi]] = varimax(AA[[mi]])$loadings
        FF[[mi]] = data %*% AA[[mi]]
        centroid = pseudoinverse(t(U) %*% U) %*% t(U) %*% FF[[mi]]
        centroid = centroid[as.integer(names(aa)), ]
    }
    else if (rotation == "promax") {
        AA[[mi]] = promax(AA[[mi]])$loadings[1:m, 1:ndim]
        FF[[mi]] = data %*% AA[[mi]]
        centroid = pseudoinverse(t(U) %*% U) %*% t(U) %*% FF[[mi]]
        centroid = centroid[as.integer(names(aa)), ]
    }
    out = list()
    mi = which.min(func)
    out$obscoord = FF[[mi]]
    rownames(out$obscoord) = rownames(data)
    out$attcoord = data.matrix(AA[[mi]])
    rownames(out$attcoord) = colnames(data)
    out$centroid = centroid
    names(cluID) = rownames(data)
    out$cluID = cluID
    out$criterion = func[mi]
    out$csize = round((table(cluID)/sum(table(cluID))) * 100, 
        digits = 1)
    out$odata = odata
    out$scale = scale
    out$center = center
    out$nstart = nstart
    class(out) = "cluspca"
    return(out)
}
            
            
    
