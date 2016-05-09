#
# Functions used to identify genetic dependencies
# by analyzing siRNA screens in cell line panels
#


library(gplots)


# Function to create a list of tables with
# common cell lines in RNAi and mutation
# data sets.

read_rnai_mutations <- function(
	rnai_file,
	func_muts_file,
	all_muts_file
	){

	rnai <- read.table(
		file=rnai_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	func_muts <- read.table(
		file=func_muts_file,
		sep="\t",
		header=TRUE,
		row.names=1
		)

	all_muts <- read.table(
		file=all_muts_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	# find the set of cell lines (ie. row names)
	# in the siRNA and mutation data sets
	common_celllines <- intersect(
		rownames(rnai),
		intersect(
			rownames(func_muts),
			rownames(all_muts)
			)
		)
	
	# return a list of tables containing only
	# the common cell lines
	return(
		list(
			rnai=rnai[common_celllines,],
			func_muts=func_muts[common_celllines,],
			all_muts=all_muts[common_celllines,]
			)
		)
}



# function to run Wilcoxon Rank Sum tests
# and Spearman Rank Correlation tests on
# each siRNA x gene combination

run_univariate_tests <- function(
	zscores,
	mutations,
	all_variants,
	alt="less"
	){
		
	zscores <- as.matrix(zscores)
	mutations <- as.matrix(mutations)
	all_variants <- as.matrix(all_variants)
	
	results <- NULL
	i <- NULL
	for(i in seq(1:length(colnames(mutations)))){
		gene <- colnames(mutations)[i]
		
		# grpA is the group of interest. e.g. mutant
		# cell lines or cell lines of a particular
		# histotype.
		grpA <- which(mutations[,i] > 0)

		# grpB includes cell lines with no reported
		# mutations at all in gene. The 'other' group
		grpB <- which(all_variants[,gene] == 0)

		j <- NULL
		for(j in seq(1:length(colnames(zscores)))){
			
			# skip if we have < 3 viability measurements
			# in one or other group
			if(min(
				length(na.omit(zscores[grpA,j])),
				length(na.omit(zscores[grpB,j]))
				) < 3){
				next
			}
						
			wilcox.p <- NA
			try(
				test <- wilcox.test(
					zscores[grpA,j],
					zscores[grpB,j],
					alternative=alt
				)
			)
			wilcox.p <- test$p.value
			
			# estimate rank serial correlation
			# as a common language effect size
			cles <- 1-(test$statistic/(
				length(zscores[grpA,j])*length(zscores[grpB,j])
				)) 
			
			marker <- colnames(mutations)[i]
			target <- colnames(zscores)[j]
			nA <- length(grpA)
			nB <- length(grpB)
			med.grpA <- median(zscores[grpA,j], na.rm=TRUE)
			med.grpB <- median(zscores[grpB,j], na.rm=TRUE)
			min.grpA <- min(zscores[grpA,j], na.rm=TRUE)
			min.grpB <- min(zscores[grpB,j], na.rm=TRUE)
			

			results <- rbind(
				results,
				c(
					marker,
					target,
					nA,
					nB,
					med.grpA,
					med.grpB,
					min.grpA,
					min.grpB,
					cles,
					wilcox.p
				)
			)
		}
	}
	
	if(is.null(nrow(results))){
		return(NULL)
	}
	
	colnames(results) <- c(
		"marker",
		"target",
		"nA",
		"nB",
		"med.grpA",
		"med.grpB",
		"min.grpA",
		"min.grpB",
		"cles",
		"wilcox.p"
	)
	
	# FDR correction of wilcoxon p-values
	wilcox_fdr <- p.adjust(
		as.numeric(results[,"wilcox.p"]),
		method="fdr"
		) 
	results <- cbind(results,wilcox_fdr)
	return(results)	
}


#
# Function to plot a heatmap of siRNA Z-scores
# stacked with a heatmap of functional mutations
# across the cell line panel
#

zscore_mutation_heatmap <- function(
	zscores,
	mutations,
	target, # siRNA target of interest
	markers, # genes to plot mutation status of
	file="z_and_mutations_plot.pdf",
	bottom_margin=20,
	right_margin=6
	){
	
	# missing values (in z-scores) are a massive pain
	# for this reason, I'm joining the rnai and mutation
	# data and na.omit()ing the lot from the start...
	
	target_mutations <- cbind(
		zscores[,target],
		mutations[,markers]
		)
	target_mutations <- na.omit(target_mutations)
	colnames(target_mutations)[1] <- target
	
	pathway_status <- apply(target_mutations[,-1],1,max)
	target_mutations <- cbind(pathway_status, target_mutations)
	colnames(target_mutations)[1] <- "pathway_status"

	
	# get the sort order for cell lines by z-score
	zscore_sort_order <- sort(
		target_mutations[,target],
		index.return=TRUE
		)
	
	# any cell lines with missing values for z-scores
	# above will have been removed from zscore_sort_order
	# need to be able to exclude these cell lines from the 
	# mutation data too.
	
	
	celllines_ordered <- rownames(target_mutations)[zscore_sort_order$ix]
	
	celllines_ordered_axis_labels <- sub(
		"_.+",
		"",
		celllines_ordered,
		perl=TRUE
		)
	
	# z-score heatmap colour palette
	breaks=seq(-4, 4, by=0.2) 
	breaks=append(breaks, 100)
	breaks=append(breaks, -100, 0)
	mycol <- colorpanel(
		n=length(breaks)-1,
		low="cyan",
		mid="black",
		high="yellow"
		)
	
	
	pdf(
		file=file,
		width=14,
		height=5
		)
	
	layout(
		matrix(c(1,2,3,3,3,3,3,3,3,3),10,1,byrow=TRUE),
		respect=FALSE
		)
	par(oma=c(0,0,0.1,0))
	
	# add z-score colour bar
	par(mar=c(0.5,9,0.5,right_margin))
	
	image(
		as.matrix(target_mutations[celllines_ordered,target]),
		xaxt="n",
		yaxt="n",
		col=mycol,
		breaks=breaks
		)
	axis(
		side=2,
		at=0.5,
		labels=paste("si", strsplit(target,"_")[[1]][1], "\nz-score", sep=""),
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)

	# add pathway status
	par(mar=c(0.5,9,0.5,right_margin))
	
	image(
		as.matrix(target_mutations[celllines_ordered,"pathway_status"]),
		xaxt="n",
		yaxt="n",
		col=c("#FFFFFF","#303030")
		)
	axis(
		side=2,
		at=0.5,
		labels="Pathway\nstatus",
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
	
	# decide bottom margin size based on number of marker genes
	#bottom_margin_size <- round(30 - (length(markers_to_print) * 1.429), digits=0)
	
	
	# add mutation status by gene
	par(mar=c(bottom_margin,9,0.5,right_margin))
	image(
		as.matrix(
			target_mutations[
				celllines_ordered,
				markers
				]
			),
		xaxt="n",
		yaxt="n",
		col=c("#FFFFFF","#A7A7A7")
		)
	# axis for cell line names
	axis(
		side=1,
		at=seq(from=0, to=1, by=1/(nrow(target_mutations)-1)),
		labels=celllines_ordered_axis_labels,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
		
	# clean up tree_gene names for printing
	markers_to_print <- NULL
	i <- NULL
	for(i in 1:length(markers)){
		markers_to_print[i] <- strsplit(markers[i],"_")[[1]][1]
	}
	
	#axis for gene names
	axis(
		side=2,
		at=seq(from=0, to=1, by=1/(length(markers)-1)),
	#	at=seq(from=0, to=1, by=1/(14-1)),
		labels=markers_to_print,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
	dev.off()
	
}




