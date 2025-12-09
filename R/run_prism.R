

#' function to validate opt.control
#' @param control a named list of parameters required to control optimization  
valid.opt.control <- function(control){
	
	ctrl <- list(maxit = 100000, maximize = FALSE, trace = 0, eps = 1e-07,
				dowarn = TRUE, tol=0, maxNA=500, n.cores=1, optimizer="MAP", sigma=2)
	namc <- names(control)
    
    if (!all(namc %in% names(ctrl)))
        stop("unknown names in opt.control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    
    if(! ctrl$optimizer %in% c("MAP","MLE"))
    		stop("unknown names of optimizer: ", ctrl$optimizer)
    
    stopifnot(length(ctrl$optimizer)==1)
    
    if(ctrl$optimizer=="MAP"){
    		if(!is.numeric(ctrl$sigma))
    			stop("sigma needs to be a numeric variable")
    		else{
    			if(ctrl$sigma<0){
    				stop("sigma needs to be positive")
    			}	
    		}	
    }
        
    return(ctrl)
}


#' function to validate gibbs.control
#' @param control a named list of parameters required to control optimization  
valid.gibbs.control <- function(control){
	ctrl <- list(chain.length=1000, burn.in=500, thinning=2, n.cores=1, seed=123, alpha=1)
	namc <- names(control)
	
	if (!all(namc %in% names(ctrl)))
        stop("unknown names in gibbs.control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control

	if(ctrl$alpha <0) stop("alpha needs to be positive")
	
	return(ctrl)
}



#' main function to run deconvolution using BayesPrism with optional chain saving
#' @param prism a prism object
#' @param n.cores number of cores to use. Default =1.
#' @param update.gibbs a logical variable to denote whether run final Gibbs sampling to update theta
#' @param save.chain logical, whether to save full Gibbs chain to HDF5. Default=FALSE.
#' @param h5.file character, path to HDF5 file for chain storage. Only used if save.chain=TRUE.
#'		If NULL and save.chain=TRUE, defaults to "gibbs_chain.h5" in current directory.
#' @param gibbs.control a list containing parameters of the Gibbs sampler
#'		chain.length: length of MCMC chain. Default=1000;
#'		burn.in: length of burn in period. Default=500;
#'		thinning: retain every # of MCMC samples after the burn in period. Default=2;
#'		n.cores: number of cores to use. Default uses n.cores in the main argument.
#'		seed: seed number for reproducibility. Default = 123.
#'		alpha: parameter of dirichlet distribution. Default=1.
#' @param opt.control a list containing parameters for the optimization step
#' @return a BayesPrism object
#' @export
run.prism <- function(prism,
					  n.cores=1,
					  update.gibbs=TRUE,
					  save.chain=FALSE,
					  h5.file=NULL,
					  gibbs.control=list(),
					  opt.control=list()
					  ){
	
	if(! "n.cores" %in% names(gibbs.control)) gibbs.control$n.cores <- n.cores
	if(! "n.cores" %in% names(opt.control)) opt.control$n.cores <- n.cores
	stopifnot(is.logical(update.gibbs) & length(update.gibbs)==1)
	stopifnot(is.logical(save.chain) & length(save.chain)==1)
	stopifnot(is.numeric(n.cores) & length(n.cores)==1)
	
	# Set default HDF5 file path if saving chains but no path provided
	if(save.chain & is.null(h5.file)) {
		h5.file <- "gibbs_chain.h5"
		cat("Chain will be saved to:", h5.file, "\n")
	}
	
	opt.control <- valid.opt.control(opt.control)
	gibbs.control <- valid.gibbs.control(gibbs.control)
	
	if(prism@phi_cellState@pseudo.min==0) 
		gibbs.control$alpha <- max(1, gibbs.control$alpha)
	
	#write mixture to disk and load each X_n to the corresponding node to save memory
	tmp.dir <- tempdir(check=TRUE)
	for(n in 1:nrow(prism@mixture)) {
		X_n <- prism@mixture[n,]
		file.name <- paste(tmp.dir, "/mixture_",n,".rdata",sep="")
		save(X_n, file= file.name)
	}
			
	#sampling cell states (cs)
	gibbsSampler.ini.cs <- new("gibbsSampler",
								reference = prism@phi_cellState,
								X = prism@mixture,
								gibbs.control = gibbs.control)

	# Run initial Gibbs with optional chain saving
	jointPost.ini.cs <- run.gibbs(gibbsSampler.ini.cs, 
								  final=FALSE,
								  save.chain=save.chain,
								  h5.file=h5.file)

	#merge over cell states to get cell type (ct) info
	jointPost.ini.ct <- mergeK(jointPost.obj = jointPost.ini.cs, 
					   		   			 map = prism@map)

	if(!update.gibbs) #only do initial Gibbs sampling
		bp.obj <- new("BayesPrism",
		 				prism = prism,
         				posterior.initial.cellState = jointPost.ini.cs,
         				posterior.initial.cellType = jointPost.ini.ct,
         				control_param = list(gibbs.control = gibbs.control, 
         									 opt.control = opt.control,
         									 update.gibbs = update.gibbs,
         									 h5.file = h5.file,
         									 save.chain = save.chain)
         			 )
	else{
		#contruct the new gibbsSampler object with updated reference
		psi <- updateReference (Z = jointPost.ini.ct@Z,
								phi_prime = prism@phi_cellType,
								map = prism@map,
								key = prism@key,
								opt.control = opt.control)
	
		gibbsSampler.update <- new("gibbsSampler",
									reference = psi,
									X = prism@mixture,
									gibbs.control = gibbs.control)
		
		# Run final Gibbs (usually doesn't save chains for this step)
		theta_f <- run.gibbs(gibbsSampler.update, 
							 final=TRUE,
							 save.chain=FALSE)
		
		bp.obj <- new("BayesPrism",
		 				prism = prism,
         				posterior.initial.cellState = jointPost.ini.cs,
         				posterior.initial.cellType = jointPost.ini.ct,
         				reference.update = psi,
         				posterior.theta_f = theta_f,
         				control_param = list(gibbs.control = gibbs.control, 
         									 opt.control = opt.control,
         									 update.gibbs = update.gibbs,
         									 h5.file = h5.file,
         									 save.chain = save.chain)
         			 )
	}
	
	unlink(tmp.dir, recursive = TRUE)
	
	cat("BayesPrism run complete.\n")
	if(save.chain & !is.null(h5.file))
		cat("Gibbs chain saved to:", h5.file, "\n")
	
	return(bp.obj) 		
}



#' function to run update gibbs sampling for BayesPrism output with initial gibbs results
#' @param bp a BayesPrism output with initial gibbs results
#' @param gibbs.control a list containing parameters of the Gibbs sampler
#'		chain.length: length of MCMC chain. Default=1000;
#'		burn.in: length of burn in period. Default=500;
#'		thinning: retain every # of MCMC samples after the burn in period. Default=2;
#'		n.cores: number of cores to use. Default =1.
#'		seed: seed number for reproducibility. Default = 123.
#'		alpha: parameter of dirichlet distribution. Default=1.
#' @param opt.control a list containing parameters for the optimization step
#' @return a BayesPrism object
#' @export
update.theta <- function(bp,
					     gibbs.control=list(),
					     opt.control=list()){
	
	#use parameters from BayesPrism run if no change
	gibbs.control.bp <- bp@control_param$gibbs.control
	opt.control.bp <- bp@control_param$opt.control

	gibbs.control.bp[names(gibbs.control)] <- gibbs.control
	opt.control.bp[names(opt.control)] <- opt.control

	opt.control <- valid.opt.control(opt.control.bp)
	gibbs.control <- valid.gibbs.control(gibbs.control.bp)

	
	#contruct the new gibbsSampler object with updated reference
	psi <- updateReference (Z = bp@posterior.initial.cellType@Z,
							phi_prime = bp@prism@phi_cellType,
							map = bp@prism@map,
							key = bp@prism@key,
							opt.control = opt.control)
	
	gibbsSampler.update <- new("gibbsSampler",
								reference = psi,
								X = bp@prism@mixture,
								gibbs.control = gibbs.control)
		
	theta_f <- run.gibbs(gibbsSampler.update, final=TRUE)
		
	bp.updated <- new("BayesPrism",
		 				prism = bp@prism,
         				posterior.initial.cellState = bp@posterior.initial.cellState,
         				posterior.initial.cellType = bp@posterior.initial.cellType,
         				reference.update = psi,
         				posterior.theta_f = theta_f,
         				control_param = list(gibbs.control = gibbs.control, 
         									 opt.control = opt.control,
         									 update.gibbs = TRUE)
         			 )
       
	return(bp.updated) 
	
}




					
