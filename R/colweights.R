

### Mandatory:
# data: any sparse or dense matrix
# weightfunction: "beta" "step" "linear" "IDF" "IDF^2"
# defualts to "IDF"

### Optional:
# binary = do you want to convert your data to binary?
# par1, par2: beta parameters or step function boundaries 
# mode: linear function max (desired density max weight, e.g. 1/k)
# convertsparse = if your matrix is dense (wide format), do you want to 
#                   convert it to a sparse matrix?
# lower = lower bound for column density trimming   (inclusive)
# upper = upper bound for column density trimming   (inclusive)


### Return:
# a Matrix-type object


### NOTE: this function assumes that multiplying by sqrt(weight) is appropriate for similarity calculation


colweights <- function (data, weightfunction='IDF',binary=True,warn=False,...) {
    par <- list(...)
    sparseinput <-is(data, 'sparseMatrix')
    require(Matrix)
    
    #####
    ##### construct Matrix object for use with Matrix package #####
    
    if (binary==T) {
        data <- data > 0
        colsum <- colSums(data)
        colprop <- colsum/nrow(data)
    }else{
        colsum <- colSums(data > 0)
        colprop <- colsum/nrow(data)
    }
 
    #####
    #####
    ##### Find density proportion of each column #####
    
    weightfunction <- as.character(weightfunction)
    
    #####
    #####
    ##### Remove columns outside your threshold (and monitor the rows) #####
    
    if( !(is.null(par$lower))) {
        
        data <- data[,which(colsum >= par$lower)]   # colsum is the sum of nonzero entries
        
        colsum <- colsum[which(colsum >= par$lower)]
        
    }
    
    if( !(is.null(par$upper))) {
        
        data <- data[,which(colsum <= par$upper)]   # colsum is the sum of nonzero entries
        
        colsum <- colsum[which(colsum <= par$upper)]
        
    }
    
    rowsum <- rowSums(data)
    
    if(min(rowsum) <= 0 & warn) {  # some rows could lose all nonzero entries when you trim columns
        #What is this for? How do I fix this?
        resp <- readline(prompt="One or more rows has zero weight. \n 
            Make sure that you fix this before continuing. \n 
            Press the ENTER key to continue. \n")
    
    }
        
        
    #####
    ##### 
    ##### Calculate column-weighted matrix & return #####
    
    if (sparseinput==F) { # if you insist on using a dense matrix
        
        if (weightfunction == "beta") {   # par1 = alpha, par2 = beta
            
            x <- seq(0,1, length=1000)
            mode.beta <- max(dbeta(x, shape1=par$par1, shape2=par$par2))
            
            colweights <- dbeta(colprop, shape1=par$par1, shape2=par$par2)
            colweights <- colweights/max(mode.beta)    # scale to (0,1) range
            colweights <- sqrt(colweights)
            return(t(t(data)/colweights))
            
        }
        
        else if (weightfunction == "step") {    # par1 = min cutoff, par2 = max cutoff
            
            return(data[,colprop > par$par1 & colprop < par$par2])
            
        }
        
        else if (weightfunction == "linear") {
            
            slope1 = 1/par$mode
            slope2 = -1/(1-par$mode)
            
            linweight <- function (density) {
                if (density < par$mode) { return(slope1*density) }
                else{return(slope2*(density-1) ) }
            }
            
            colweights <- sapply(colprop, linweight)
            colweights <- sqrt(colweights)
            return(t(t(data)/colweights))
            
        }
        
        else if (weightfunction == "IDF") {
            
            # IDF column weighting = log( N/ 1+density )
            data.idf <- log(nrow(data)/(1 + colsum))
            data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)
            
            # multiply each column by its IDF weight
            data.tfidf <- crossprod(t(data), data.idf.diag)
            return(data.tfidf)
            
            # Row normalize
            # data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
            # data.tfidf.rn <- data.tfidf/ rowSums(data.tfidf)
            # return(data.tfidf.rn)
            
        }
        
        else if (weightfunction == "IDF^2") {
            
            # IDF column weighting = log( N/ 1+density )
            data.idf <- (log(nrow(data)/(1 + colsum)))^2
            
            # Multiply each column by its IDF weight
            data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)
            data.tfidf <- crossprod(t(data), data.idf.diag)
            return(data.tfidf)
            
            # Row normalize
            # data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
            # data.tfidf.rn <- data.tfidf/ rowSums(data.tfidf)
            # return(data.tfidf.rn)
            
        }
        
        else if (weightfunction == "none") {
            
            return(data)
            
        }
        
        else {stop("Pick a valid weight method.")}
        
    }
    
    else { # sparse matrix calculations
        
        if (weightfunction == "beta") {   # par1 = alpha, par2 = beta
            
            x <- seq(0,1, length=1000)
            mode.beta <- max(dbeta(x, shape1=par$par1, shape2=par$par2))
            
            colweights <- dbeta(colprop, shape1=par$par1, shape2=par$par2)
            colweights <- colweights/max(mode.beta)    # scale to (0,1) range
            colweights <- sqrt(colweights)
            return(t(t(data)/colweights))
            
        }
        
        else if (weightfunction == "step") {    # par1 = min cutoff, par2 = max cutoff
            
            return(data[,colprop > par$par1 & colprop < par$par2])
            
        }
        
        else if (weightfunction == "linear") {
            
            slope1 = 1/par$mode
            slope2 = -1/(1-par$mode)
            
            linweight <- function (density) {
                if (density < par$mode) { return(slope1*density) }
                else{return(slope2*(density-1) ) }
            }
            
            colweights <- sapply(colprop, linweight)
            colweights <- sqrt(colweights)
            return(t(t(data)/colweights))
            
        }
        
        else if (weightfunction == "IDF") {
            
            # IDF column weighting = log( N/ density )
            data.idf <- log(nrow(data)/(colsum))
            
            # Multiply each column by its IDF weight
            data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)
            data.tfidf <- crossprod(t(data), data.idf.diag)
            return(data.tfidf)
            
            # Row normalize
            # data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
            # data.tfidf.rn <- data.tfidf/ rowSums(data.tfidf)
            # return(data.tfidf.rn)
            
        }
        
        else if (weightfunction == "IDF^2") {
            
            # IDF column weighting = log( N/ density )
            data.idf <- (log(nrow(data)/(colsum)))^2
            
            # Multiply each column by its IDF weight
            data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)
            data.tfidf <- crossprod(t(data), data.idf.diag)
            return(data.tfidf)
            
            # Row normalize
            # data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
            # data.tfidf.rn <- data.tfidf/ rowSums(data.tfidf)
            # return(data.tfidf.rn)
            
        }
        
        else if (weightfunction == "none") {
            
            return(data)
            
        }
        
        else {stop("Pick a valid weight method.")}
        
    }
    
}
