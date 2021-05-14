source("Scripts/Experiment.Core.R")

# Load database files -----------------------------------------------------

# load files from a path i.e ( Data )
# a list of dataframes will be enlisted in a variable ( dataset )

files <- list.files(path = "Data",pattern = ".txt", full.names = TRUE)
dataset <- lapply(files, vroom, delim = "\t") 
filenames <- list.files(path = "Data",
                        pattern = ".txt") %>% gsub(".txt.*", "", .)
names(dataset) <- filenames


# Prepare assay for Experiment.Core ---------------------------------------

# Merge all dataframes enlisted in the variable
# This creates a NAMED count matrix from the dataframes

PDAC.counts <- Reduce(merge,dataset) 
rownames(PDAC.counts) <- PDAC.counts$name
PDAC.counts <- PDAC.counts %>% dplyr::select(-name) %>% as.matrix

# Metadata is a dataframe supplementary to the count matrix
# It consists of SAMPLES their possible groupings and any possible supporting information
# Current design of Experiment.Core only supports experiment of this setup 

PDAC.meta <- colnames(PDAC.counts) %>% as.data.frame
colnames(PDAC.meta) <- "SAMPLES"
PDAC.meta <- PDAC.meta %>%  mutate(Group = as.factor(casefold(gsub(".*?([A-Z,a-z]+).*", "\\1", SAMPLES)
                                                              ,upper = TRUE))) %>%
  mutate(FIS.CON = as.factor(paste0("CON",if_else(.$Group == 'FIS',gsub(".*?([0-9]+).*", "\\1", SAMPLES), '0'))))



# Create an object for Experiment.Core ------------------------------------

# Listing the prior raw data prepared for Experiment.Core to make it easier to remove if needed

tmp <- ls()

# Instantiates an R6 Object [ExperimentCore] , a custom class aimed to treat experimental data
# with one or more observations ( i.e. Control and treatment or here Fis or NC)
# Mostly it is an encapsulation of custom functions and functions from some packages; consolidated
# to provide a more straightforwared pipeline

# CreateDGE [bool] - for setting if the instantiation creates the DGEList from start
# meta [dataframe] -  the defined metadata  which is supplementary to the 
# assay [matrix] - a named count matrix (ENSMBL for now)
# controls [dataframe] - a reference to the observed group from the metadata
# Complete structure of ExperimentCore can be seen through : str(ExperimentCore)

PDAC.Experiment <- ExperimentCore$new(createDGE = TRUE, meta = PDAC.meta, 
                                      assay = PDAC.counts, controls = PDAC.meta$Group)




# Starts the analysis pipeline --------------------------------------------

# Creates a design matrix needed for Differential Expression
# dependency : edgeR

# remove.lowcounts [bool] - enabling removal of low counts throught lib.size observation ( with the help 
# of edgeR )
# custom.design [model.matrix] - a parameter when one wants to provide a custom design matrix - may need to modify
# the metadata [this parameter is optional]

PDAC.Experiment$designExperiment(remove.lowcounts = TRUE)

# log2 [bool] - to perform a log2 transformed data default is TRUE  
# of edgeR )
# lib.size [bool] - to consider lib.size during normalization default is TRUE
# prior.count [integer] - to add the counts during normalization default is 1
# method [character] - selected normalization method default is TMM (options: TMM",
# "RLE","upperquartile")
# if parametes are left unassigned it will assume the default values

PDAC.Experiment$normalizeData()

# Performs dispersion analysis - NB and Quasi-likelihood , parameters added to $dge object under
# instantiated ExperimentCore (i.e PDAC.Exmperiment$dge)
# plot will be printed along 

PDAC.Experiment$estimateDisperion()

# Performs differential expression analysis between 2 observations
# variable1 and variable 2 [character]
# threshold [numeric] - the limit at which genes will be classified as DE or not
# stores result in $test - a list of all test done (QL , LRT)
# test can be accessed through [name_of_ExperimentCore]$test[['LRT]] 
# i.e. PDAC.Experiment$test[['LRT']]

PDAC.Experiment$contrastVariables(variable1 = 'FIS', variable2 = 'NC')
PDAC.Experiment$performDifferential(threshold = 1.5)

# Originally an MA plot will be created but if user tend to prefer a ggplot based volcano plot
# the user can prompt the command below
# accepts the following parameters:
# test [character] - the either QL or LRT ( default : LRT) , the later has consideration of the threshold
# p.threshold [numeric] - needed for setting a marker line for p.adjusted value to determine DE genes 
# default is 0.05 ( genes exceeding the p.threshold)
# logfc.threshold [numeric] - needed for setting a marker line for log2fc value to determine DE genes 
# default is 1.5 ( genes exceeding the logfc.threshold)

# plot will be saved in $plots - a list of ExperimentCore plots, can be accessed by name or index
# ie. PDAC.Experiment$plots[1]
# names(PDAC.Experiment$plots) will show all the names of plots
# invoke PDAC.Experiment$plots if you want to see all plots

PDAC.Experiment$volcanoPlot()

# Example:
PDAC.Experiment$plots

# Or
PDAC.Experiment$plots[[1]]

# Or (as observed from the plot list)
PDAC.Experiment$plots[['VolcanoPlot1']]

# Performs aggregate statistics (MEAN,MEDIAN, QUANTILE of prob q, SD, VAR )
# q [integer] - represents the probability in % i.e 10 (decile), 25 (quartile) etc
# result will be stored in object$stat ie PDAC.Experiment$stat
PDAC.Experiment$calcStatistics(q = 20)

# Creates a ranked list of genes based from their selected statistic
# stat [character] - the selected stastictic to rank default is logFC (others can be 
# observed from $test$table i.e. logFC, logCPM, PValue )
# result will be stored in object$ranked.genes ie PDAC.Experiment$ranked.genes

PDAC.Experiment$createRankedList()

# Creates an dot plot and barplot of normalized count of DE genes showing contrast between
# the observations
# n.genes [integer] -  number of DE genes to be shown
# result will be stored in object$plots ie PDAC.Experiment$plots[[index]]

PDAC.Experiment$plotExpression(n.genes = 10)

#names(PDAC.Experiment$plots)
PDAC.Experiment$plots[['Expression2']]
PDAC.Experiment$plots[['Expression3']]


# Initialize Gene Ontology database by defining an initial GO enrichment with given terms
# terms [charcter] or vector of [characters] whose values range from BP, MF and CC
# terms can be given as c('BP','MF', 'CC')
# database will be stored in object$GO.Database
# enrichment in obejct$GO.Analysis, obejct$GO.Enrichment and object$Result 
# GO plots will be added on

PDAC.Experiment$initializeGO(terms = 'BP')

# n [numeric] - the top n of the GO terms default is 15 
PDAC.Experiment$enrichGO(GO.Term = 'BP')

#names(PDAC.Experiment$plots)
# NOTE : GoPlot shows the
PDAC.Experiment$plots[['GoPlot4']]
PDAC.Experiment$plots[['Enrichment5']]


# Add plots for another selected terms/pathways
# The listed pathways/terms are in object$GO.Analysis
PDAC.Experiment$addGOPlot(GO.Pathway = 'single strand break repair')


# Adding GO Analysis based from another term (options : MF, BP, CC)
# This is to run additional GO Analysis on terms not included in prior analysis
PDAC.Experiment$runGO(terms = 'MF')
PDAC.Experiment$enrichGO(GO.Term = 'MF')

# initialize GSEA database 
# pathway.db [character] - the directory of gmt file to be studied
# name [character] - name to be assigned to the database on its listing to object$GSEA.Database
PDAC.Experiment$initializeGSEA()

# name [character] - the name database to be used for GSEA
# perm [integer] - the number of permutations to be used in p val calculation default is 10000
PDAC.Experiment$performGSEA()

# name [character] - the name database to be used for Gene Enrichment
# n [integer] - the top n of the GSEA pathways default is 15 
PDAC.Experiment$plotEnrichedPathways()

# Creates a GSEA enrichment plot
# database [character] -  name of database to be used in the enrichment
# pathway [character] - name of the desired pathway
# perm [integer] - the number of permutations to be used in p val calculation default is 10000
PDAC.Experiment$plotGSEA()
PDAC.Experiment$plotGSEA(pathway = 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE')
PDAC.Experiment$plotGSEA(pathway = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB')

PDAC.Experiment$plots[['Enrichment9']]

for(i in names(PDAC.Experiment$plots)){
  print(PDAC.Experiment$plots[[i]])
  ggsave(paste0("Plots/",i,".png"))
  
}
