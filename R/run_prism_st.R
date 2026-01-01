

run.prism.st <- function(prism,
					  	 n.cores=1,
					  	 save.chain="none",
					  	 h5.file=NULL,
					  	 compute.z.cv=FALSE,
					  	 gibbs.control=list(),
					  	 opt.control=list(optimizer="MLE")){
	
	if(! "n.cores" %in% names(gibbs.control)) gibbs.control$n.cores <- n.cores
	if(! "n.cores" %in% names(opt.control)) opt.control$n.cores <- n.cores
	stopifnot(is.numeric(n.cores) & length(n.cores)==1)
	stopifnot(is.character(save.chain) & length(save.chain)==1)
	stopifnot(save.chain %in% c("none", "theta", "all"))
	stopifnot(is.logical(compute.z.cv) & length(compute.z.cv)==1)
	
	# Set default HDF5 file path if saving chains but no path provided
	if(save.chain != "none" & is.null(h5.file)) {
		h5.file <- "gibbs_chain.h5"
		cat("Chain will be saved to:", h5.file, "\n")
	}
	
	opt.control <- valid.opt.control(opt.control)
	gibbs.control <- valid.gibbs.control(gibbs.control)
			
	#sampling cell states (cs)	
	gibbsSampler.ini <- new("gibbsSampler",
								reference = prism@phi_cellState,
								X = prism@mixture,
								gibbs.control = gibbs.control)

	jointPost.ini <- run.gibbs(gibbsSampler.ini, final=FALSE, 
							   save.chain=save.chain, h5.file=h5.file, compute.z.cv=compute.z.cv)

	#perform MLE estimate for gamma
	psi <- updateReference (Z = jointPost.ini@Z,
							phi_prime = prism@phi_cellState,
							map = prism@map,
							key = prism@key,
							opt.control = opt.control)
	
	gibbsSampler.update <- new("gibbsSampler",
								reference = psi,
								X = prism@mixture,
								gibbs.control = gibbs.control)
	
	jointPost.update <- run.gibbs(gibbsSampler.update, final=FALSE,
								  save.chain="none", compute.z.cv=compute.z.cv)

	#merge over cell states to get cell type (ct) info
	jointPost.update.ct <- mergeK(jointPost.obj = jointPost.update, 
					   		   	  map = prism@map)
	
	bp.obj <- new("BayesPrismST",
		 		   prism = prism,
		 		   posterior.cellState = jointPost.update,
		 		   posterior.cellType = jointPost.update.ct,
		 		   reference.update = psi,
		 		   control_param = list(gibbs.control = gibbs.control, 
         								opt.control = opt.control,
         								save.chain = save.chain,
         								h5.file = h5.file))
	
	cat("BayesPrism.st run complete.\n")
	if(save.chain != "none" & !is.null(h5.file))
		cat("Gibbs chain saved to:", h5.file, "\n")
         			 
	return(bp.obj) 		
}






