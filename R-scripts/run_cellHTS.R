#
# Script to process a single siRNA screen
# using cellHTS2
#

# You need to change this to point to the
# actual directory the data are in
setwd("~/Dropbox/DrJCampbell/identifying-genetic-dependencies/cellHTSanalysis")

# load the cellHTS2 package
require(cellHTS2)

# read the plate list file defining the filenames
# corresponding to each screen plate and replicate. This
# creates a cellHTS object we will call x (see Note 1)
x <- readPlateList(
	filename="platelist_p3r3.txt",
	name="CGDsExample",
	path="./"
	)

# add information to the cellHTS object to indicate the
# positions of the controls on the plates, the screen 
# description file and the screen log file
x <- configure(
	x,
	descripFile="description.txt",
	confFile="kinome_384_plateconfig.txt",
	path="./"
	)

# add information to the cellHTS2 object to indicate
# which genes are targeted by the siRNAs in each well
x <- annotate(
	x,
	geneIDFile="kinome_library.txt",
	path="./"
	)

# perform normalisation of each plate in the screen. Log
# transform and center data to the median intensity of
# each plate
xn <- normalizePlates(
	x,
	scale="multiplicative",
	log=TRUE,
	method="median",
	varianceAdjust = "none",
	posControls="siplk1",
	negControls="sicon1|sicon2|allstar",
	)

# Scale the well intensities to the median absolute
# deviation across the library to produce Z-scores
xsc <- scoreReplicates(
	xn,
	method="zscore",
	sign="+"
	)

# Summarise the Z-scores as medians
xsc <- summarizeReplicates(
	xsc,
	summary="median"
	)

# Write the cellHTS top table providing a summary of the
# screen
summary_info <- getTopTable(
	list(
		"raw"=x,
		"normalized"=xn,
		"scored"=xsc
		),
	file="TopTable.txt"
	)

# configure the HTML report created by writeReport()
setSettings(
	list(
		plateList=list(
			reproducibility=list(
				include=TRUE,
				map=TRUE
				),
			intensities=list(
				include=TRUE,
				map=TRUE)
				),
		screenSummary=list(
			scores=list(
				range=c(-20, 10),
				map=TRUE
				)
			)
		)
	)

# Produce a detailed HTML report in a subdirectory
writeReport(
	raw=x,
	normalized=xn,
	scored=xsc,
	outdir="./report",
	force=TRUE,
	posControls="siplk1",
	negControls="sicon1|sicon2|allstar",
	mainScriptFile="../R-scripts/run_cellHTS.R"
	)


# For convenience, write Z-scores to a file that can be
# joined with further screens for analysis 
plates <-plate(xsc)
wells <- well(xsc)
scores <- Data(xsc)
compounds <- geneAnno(xsc)
combinedz <- data.frame(
	compound=compounds,
	plate=plates,
	well=wells,
	zscore=scores
	)

write.table(
	combinedz,
	"zscores.txt",
	sep="\t",
	quote=FALSE,
	row.names=FALSE
	)

