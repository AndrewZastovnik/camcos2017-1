#' Checking cluster accuracy
#'
#' In the results table, true clusters are the rows, guessed clusters are the
#' columns
#'
#' @param indices: Labels from clustering
#' @param trueLabels: Ground truth labels.
#'
#' @return List of 3. Confusion matrix, label mapping, accuracy.
#' @export
clustercheck<-function(indices,trueLabels){
  # find the true cluster sizes
  K <- length(unique(trueLabels))
  planeSizes <-rep(0,K)
  k <- 1
  for(m in unique(trueLabels)){
    planeSizes[k] = sum(trueLabels == m)
    k = k+1
  }
  # num is a little confusing... matrix!
  num <- matrix(0,nrow=K,ncol=length(unique(indices)))
  k<- 1
  for(n in unique(trueLabels)){
    j=1
    for(m in sort(unique(indices))){
      # fill the confusion matrix
      num[k,j] = sum((indices==m)*(trueLabels==n))
      j=j+1
    }
    k = k+1
  }
  
  results = maximum_number_of_correctly_classified_points(num)
  p = results[[1]]/sum(planeSizes)
  return(list(confusion_matrix=num,
              accuracy = p, 
              assignments= results[[2]]))
}

maximum_number_of_correctly_classified_points <-  function(num){
  # this function finds optimal pairing of clusters and classes
  n <- nrow(num)
  m<- ncol(num)
  index_row <- rep(FALSE,n)
  index_column <- rep(FALSE,m)
  # assignments of -1 mean that cluster was given no class
  assignment <- rep(-1,n)
  count <- 0
  while(count <= n){
    # make a temp confusion matrix the we will remove obviously optimal
    # combinations of rows and columns
    num_temp = as.matrix(num[index_row==0,index_column==0])
    if (sum(num_temp)==0){
      # end the the loop if there are no nonzero entries left in the 
      #confusion matrix
      break
    }
    # the rows of the original confusion matrix that haven't been 
    # assigned yet
    rows <- (1:n)[index_row ==0]
    # the columns that haven't been assigned
    columns <- (1:m)[index_column == 0]
    # entries of our confusion matrix that are the max in their
    # row and column winners[[1]] is rows winners[[2]] are the columns
    winners = find_clear_winners(num_temp,rows,columns)
    assignment[winners[[1]]] = winners[[2]]
    # make sure we remove the winners from the matrix next time
    index_column <- index_column | (1:m %in% winners[[2]])
    index_row <- index_row |  (1:n %in% winners[[1]])
    if(length(columns) ==0){
      # end if there are no more clusters to assign
      # should be handeled by the other check but whatevere
      break
    }
    count <- count + 1
  }
  n <- sum(diag(num[(1:n)[assignment!=-1],assignment[assignment!=-1]]))
  return(list(assignment,n))
}

find_clear_winners<- function(num,rows,columns){
  # finds the row and column of enties that are 
  # the largest entries in the 
  max_row <- max.col(num)
  max_col <- max.col(t(num[, max_row]))
  index <- max_col == 1:nrow(num)
  row_winners <-rows[(1:nrow(num))[index]]
  column_winners <- columns[max_row[index]]
  return(list(row_winners,column_winners))
}
