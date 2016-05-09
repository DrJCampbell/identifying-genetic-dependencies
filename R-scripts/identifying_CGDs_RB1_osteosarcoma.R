#
#
#

# Update the line below to point to where you cloned
# the git repository
setwd("~/Dropbox/DrJCampbell/identifying-genetic-dependencies")

source("./R-scripts/identifying_CGDs_library.R")

sirna_screens_file <- "./siRNA-data/Osteosarcoma_kinome_screens.txt"
rb_pathway_func_muts_file <- "./mutation-data/combined_exome_cnv_func_muts_RBpathway_160509.txt"
rb_pathway_all_muts_file <- "./mutation-data/combined_exome_cnv_all_muts_RBpathway_160509.txt"

# read in and combine the siRNA and mutation
# data files
kinome_rb_muts <- read_rnai_mutations(
	rnai_file=sirna_screens_file,
	func_muts_file=rb_pathway_func_muts_file,
	all_muts_file=rb_pathway_all_muts_file
	)

# run the association tests
kinome_rb_mut_associations <- run_univariate_tests(
	zscores=kinome_rb_muts$rnai,
	mutations=kinome_rb_muts$func_muts,
	all_variants=kinome_rb_muts$all_muts,
	alt="less"
	)

# write out the associations table
write.table(
	kinome_rb_mut_associations,
	"./results/kinome_rb_mut_associations.txt",
	sep="\t",
	col.names=TRUE,
	row.names=FALSE,
	quote=FALSE
	)


zscore_mutation_heatmap(
	zscores=kinome_rb_muts$rnai,
	mutations=kinome_rb_muts$func_muts,
	target="DYRK1A_ENSG00000157540",
	markers=c(
		"CCND1_595_ENSG00000110092",
		"CDKN2A_1029_ENSG00000147889",
		"RB1_5925_ENSG00000139687"
		),
	file="./results/RB_pathway_mutations_and_DYRK1A_Zscores.pdf",
	bottom_margin=20,
	right_margin=6
	)
