#all operations for the jointPost S4 class


#' constructor of the class jointPost
#'
#' @param mixtureID, a character vector to represent identifiyers of bulk samples
#' @param cellType, a character vector to represent identifiyers of cell types
#' @param geneID, a character vector to represent identifiyers of genes
#' @param gibbs.list, list containing the posterior means of gibbs sampling, with each element containing Z_n and theta_n 
#'
#' @return a jointPost S4 object
newJointPost <- function(bulkID,
						 geneID,
						 cellType,
						 gibbs.list){
		
	N <- length(bulkID)
	G <- length(geneID)
	K <- length(cellType)
	
	stopifnot(length(gibbs.list) == N)
						 	
	Z <- array(NA, 
			   dim = c(N,G,K),
			   dimnames=list(bulkID, geneID, cellType))

	# store CV for Z (coefficient of variation) — ensure present with same dims
	Z.cv <- array(NA,
			   dim = c(N,G,K),
			   dimnames=list(bulkID, geneID, cellType))
	
	theta <- matrix(NA, 
					nrow = N, ncol= K, 
					dimnames = list(bulkID, cellType))
	
	theta.cv <- matrix(NA, 
					   nrow = N, ncol= K, 
					   dimnames = list(bulkID, cellType))
	
	for (n in 1:N) {
		Z[n,,] <- gibbs.list[[n]]$Z_n
		theta[n,] <- gibbs.list[[n]]$theta_n
	}

	# populate Z.cv if available in gibbs.list entries, otherwise leave NA
	# Z.cv_n will be NULL if compute.cv=FALSE
	if(!is.null(gibbs.list[[1]]$Z.cv_n)){
		for (n in 1:N) Z.cv[n,,] <- gibbs.list[[n]]$Z.cv_n
	}

	# theta.cv is always computed
	if(!is.null(gibbs.list[[1]]$theta.cv_n)){
		for (n in 1:N) theta.cv[n,] <- gibbs.list[[n]]$theta.cv_n
	}
	
	constant <- sum(unlist(lapply(gibbs.list, '[[', "gibbs.constant")))
	 
	new("jointPost", Z = Z, Z.cv = Z.cv, theta = theta, theta.cv = theta.cv, constant = constant)
}


#' constructor of the class jointPost
#'
#' @param mixtureID, a character vector to represent identifiyers of bulk samples
#' @param cellType, a character vector to represent identifiyers of cell types
#' @param geneID, a character vector to represent identifiyers of genes
#' @param gibbs.list, list containing the posterior means of gibbs sampling, with each element containing Z_n and theta_n 
#'
#' @return a jointPost S4 object
newThetaPost <- function(bulkID,
						 cellType,
						 gibbs.list){
		
	N <- length(bulkID)
	K <- length(cellType)
	
	stopifnot(length(gibbs.list) == N)
						 		
	theta <- matrix(NA, 
					nrow = N, ncol= K, 
					dimnames = list(bulkID, cellType))
	
	theta.cv <- matrix(NA, 
					   nrow = N, ncol= K, 
					   dimnames = list(bulkID, cellType))
	
	for (n in 1:N) {
		theta[n,] <- gibbs.list[[n]]$theta_n
		theta.cv[n,] <- gibbs.list[[n]]$theta.cv_n
	}
	
	new("thetaPost", theta = theta, theta.cv = theta.cv)
}


#' function to marginalize K, e.g. cell states within each cell type
#' 
#' @param jointPost.obj a jointPost object
#' @param map a list of the format list(cell.type1=c(cell.stateA, cell.stateB), ...)
#' @return a new jointPost object 
mergeK <- function(jointPost.obj,
				   map){
								
	bulkID <- dimnames(jointPost.obj@Z)[[1]]
	geneID <- dimnames(jointPost.obj@Z)[[2]]
	cellType <- dimnames(jointPost.obj@Z)[[3]]
	cellType.merged <- names(map)
	
	N <- length(bulkID)
	G <- length(geneID)
	K <- length(cellType)
	K_merged <- length(cellType.merged)

	stopifnot(length(unlist(map)) == K)
					 	
	Z <- array(NA, 
			   dim = c(N, G, K_merged),
			   dimnames=list(bulkID, geneID, cellType.merged))
	
	theta <- matrix(NA, 
					nrow = N, ncol= K_merged, 
					dimnames = list(bulkID, cellType.merged))

	# prepare Z.cv for merged results (same dims as Z)
	Z.cv <- array(NA,
			   dim = c(N, G, K_merged),
			   dimnames=list(bulkID, geneID, cellType.merged))
	
	theta.cv <- matrix(NA, 
					   nrow = N, ncol= K_merged, 
					   dimnames = list(bulkID, cellType.merged))
	
	#merge across cellType.merged
	for(k in 1:K_merged){
		cellType.merged.k <- names(map)[k]
		cellTypes.k <- map[[k]]
		if(length(cellTypes.k)==1) {
			#skipping summation. assign value directly
			Z[,,cellType.merged.k] <- jointPost.obj@Z[,,cellTypes.k, drop=F]
			theta[,cellType.merged.k] <- jointPost.obj@theta[,cellTypes.k, drop=F]
			# copy CVs directly if present
			if(!is.null(dim(jointPost.obj@Z.cv)))
				Z.cv[,,cellType.merged.k] <- jointPost.obj@Z.cv[,,cellTypes.k, drop=F]
			if(!is.null(dim(jointPost.obj@theta.cv)))
				theta.cv[,cellType.merged.k] <- jointPost.obj@theta.cv[,cellTypes.k, drop=F]
		}
		else {
			#marginalize cellTypes.k
			Z[,, cellType.merged.k] <- rowSums(jointPost.obj@Z[,,cellTypes.k, drop=F], dims=2)
			theta[,cellType.merged.k] <- rowSums(jointPost.obj@theta[,cellTypes.k, drop=F])
			
			# merge CVs for Z: convert CV->variance, sum variances (assume independence), then convert back to CV
			if(!is.null(dim(jointPost.obj@Z.cv))){
				# jointPost.obj@Z is N x G x K; get the subset for the cellTypes.k
				Z_means <- jointPost.obj@Z[,,cellTypes.k, drop=F] # N x G x ksub
				Z_cvs <- jointPost.obj@Z.cv[,,cellTypes.k, drop=F]
				# compute variances = (cv * mean)^2, sum across ksub
				var_sum <- apply((Z_cvs * Z_means)^2, c(1,2), sum)
				mean_sum <- Z[,, cellType.merged.k]
				# avoid divide by zero
				cv_merged <- array(NA, dim = dim(mean_sum))
				zero_mean <- mean_sum == 0
				cv_merged[!zero_mean] <- sqrt(var_sum[!zero_mean]) / mean_sum[!zero_mean]
				cv_merged[zero_mean] <- NA
				Z.cv[,,cellType.merged.k] <- cv_merged
			}
			
			# merge CVs for theta: same approach
			if(!is.null(dim(jointPost.obj@theta.cv))){
				theta_means <- jointPost.obj@theta[,cellTypes.k, drop=F] # N x ksub
				theta_cvs <- jointPost.obj@theta.cv[,cellTypes.k, drop=F]
				# compute variances = (cv * mean)^2, sum across ksub
				var_sum <- rowSums((theta_cvs * theta_means)^2)
				mean_sum <- theta[, cellType.merged.k]
				# avoid divide by zero
				cv_merged <- rep(NA, length(mean_sum))
				zero_mean <- mean_sum == 0
				cv_merged[!zero_mean] <- sqrt(var_sum[!zero_mean]) / mean_sum[!zero_mean]
				cv_merged[zero_mean] <- NA
				theta.cv[,cellType.merged.k] <- cv_merged
			}
		}
	}
	 
	new("jointPost", Z = Z, Z.cv = Z.cv, theta = theta, theta.cv = theta.cv)
	
}