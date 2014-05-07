rescale.univariate.model = function(model, reference.model) {

	## Intercept user input
	if (check.univariate.model(model)!=0) stop("argument 'model' expects a univariate model")
	if (check.univariate.model(reference.model)!=0) stop("argument 'reference.model' expects a univariate model")

	## Get variables for rescaling
	binsize.old = diff(model$coordinates$start)[1]
	numbins.old = length(model$reads)
	binsize.new = diff(reference.model$coordinates$start)[1]
	
	# Check if rescaling is possible
	if (binsize.old %% binsize.new !=0 ) {
		stop("Cannot rescale model to binsize ",binsize.new,". The old binsize ",binsize.old," has to be a multiple of the new binsize.")
	}
	rescalef = binsize.old / binsize.new

	## Prepare output model
	new.model = model
	new.model$coordinates = NULL
	new.model$reads = NULL
	new.model$posteriors = NULL
	
	## Get states
	states = get.states(model)

	## Go through each chromosome in the input model
	for (chr in levels(model$coordinates$chrom)) {
		mask = model$coordinates$chrom == chr
		ref.mask = reference.model$coordinates$chrom == chr
		numbins.new = length(which(ref.mask))
		# Coordinates
		new.model$coordinates = rbind(new.model$coordinates, reference.model$coordinates[ref.mask,])
		# Rescale posteriors
		post.new = matrix(rep(model$posteriors[mask,], each=rescalef), ncol=3)[1:numbins.new,]
		new.model$posteriors = rbind(new.model$posteriors, post.new)
# 		# Rescale reads by repeating
# 		reads.new = rep(model$reads[mask], each=rescalef)[1:numbins.new]
		# Rescale reads by drawing from distribution
		reads.new = NULL
		rle = rle(states[mask])
		for (i1 in 1:length(rle$values)) {
			if (rle$values[i1]==0) {
				reads.new = c(reads.new, rzinbinom(rle$lengths[i1]*rescalef, w=model$weights[1], size=model$distributions[2,'r'], prob=model$distributions[2,'p']))
			} else {
				reads.new = c(reads.new, rnbinom(rle$lengths[i1]*rescalef, size=model$distributions[3,'r'], prob=model$distributions[3,'p']))
			}
		}
		new.model$reads = c(new.model$reads, reads.new[1:numbins.new])
	}

	return(new.model)
}
