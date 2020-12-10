# Required R packages
library(ape)
library(phytools)
# load treespace and packages for plotting:
library(treespace)
library(phylogram)
library(phangorn)
library(seqinr)
library(adegraphics)
library(adegenet)
library(apTreeshape)
library(ggtree)

# Set seed for reproducibility 
set.seed(1)

# Load metadata file
#SRA_metadata <- read.csv("Salmonella_outbreak_SRA_metadata.csv", header = FALSE)

# Read in phylogenetic trees
lyve_tree <- read.tree(file = "Pipeline_results/Old_results/Ecoli_romaine_outbreak/exported_trees/lyveset.newick")


# kSNP3 tree.NJ.tre, tree.ML.tre, tree.core.tre, tree.parsimony.tre
ksnp_tree <- read.tree(file = "Pipeline_results/Old_results/Ecoli_romaine_outbreak/exported_trees/ksnp3.newick")

# Cfsan
cfsan_tree <- read.tree(file = "Pipeline_results/Old_results/Ecoli_romaine_outbreak/exported_trees/cfsan.newick")

# Enterobase
#enterobase_tree <- read.tree(file = "Pipeline_results/Old_resultsE_coli_romaine_outbreak/etoki/enterobase_SRA_phylo_tree.nwk")
enterobase_tree <- read.tree(file = "Pipeline_results/Old_results/Ecoli_romaine_outbreak/exported_trees/enterobase.newick")


# Combine trees
combined_trees <- c(lyve_tree,ksnp_tree,cfsan_tree,enterobase_tree)
# Combine trees from single dataset into vector
dataset1_tree_vector <- c(lyve_tree,ksnp_tree,cfsan_tree,enterobase_tree)
dataset1_tree_vector <- c(as.phylo(lyve_tree),as.phylo(ksnp_tree),as.phylo(cfsan_tree),as.phylo(enterobase_tree))


# From this point on, you have all of the phylogenetic trees loaded 
#   and can perform any further analysis you wish. 
# Alot of the code below is still experimental and needs improvement.

#
##
### Code for subsetting trees with unmatched nodes
## Need to automate this 
# 

## Check for sample matches, between each tree
all_SRA_to_drop = c()

SRA_to_drop <- unique(enterobase_tree$tip.label[! enterobase_tree$tip.label %in% cfsan_tree$tip.label])
all_SRA_to_drop = c(all_SRA_to_drop,SRA_to_drop)
SRA_to_drop <- unique(enterobase_tree$tip.label[! enterobase_tree$tip.label %in% ksnp_tree$tip.label])
all_SRA_to_drop = c(all_SRA_to_drop,SRA_to_drop)
SRA_to_drop <- unique(cfsan_tree$tip.label[! cfsan_tree$tip.label %in% enterobase_tree$tip.label])
all_SRA_to_drop = c(all_SRA_to_drop,SRA_to_drop)

SRA_to_drop <- unique(lyve_tree$tip.label[! lyve_tree$tip.label %in% ksnp_tree$tip.label])
SRA_to_drop <- unique(ksnp_tree$tip.label[! ksnp_tree$tip.label %in% lyve_tree$tip.label])

all_SRA_to_drop <- unique(all_SRA_to_drop)

lyve_tree <- drop.tip(combined_trees[[1]], all_SRA_to_drop)
ksnp_tree <- drop.tip(combined_trees[[2]], all_SRA_to_drop)
cfsan_tree <- drop.tip(combined_trees[[3]], all_SRA_to_drop)
enterobase_tree <- drop.tip(combined_trees[[4]], all_SRA_to_drop)

lyve_tree_rooted <- root(lyve_tree,1, r = TRUE)
ksnp_tree_rooted <- root(ksnp_tree,1, r = TRUE)
cfsan_tree_rooted <- root(cfsan_tree,1, r = TRUE)
enterobase_tree_rooted <- root(enterobase_tree,1, r = TRUE)


#
##
###
#### Info on tree
###
##
#
# Quick summary of #nodes, tips, branch lengths, etc
summary(lyve_tree)
sum(lyve_tree$edge.length)

summary(cfsan_tree)
sum(cfsan_tree$edge.length)

summary(enterobase_tree)
sum(enterobase_tree$edge.length)

summary(ksnp_tree)
sum(ksnp_tree$edge.length)

#
### Calculate co-speciation (RF distance) between trees
#

## Need to:
# Make function to automate these pairwise comparisons within a vector of trees
# Save results in table
# Compare trees with cospeciation

# Test using Robinson-Foulds metric (RF)
cospeciation(ksnp_tree, cfsan_tree, distance = c("RF"), method=c("permutation"), nsim = 1000)
cospeciation(ksnp_tree, enterobase_tree, distance = c("RF"), method=c("permutation"), nsim = 1000)
cospeciation(lyve_tree,ksnp_tree, distance = c("RF"), method=c("permutation"), nsim = 1000)
cospeciation(lyve_tree,cfsan_tree, distance = c("RF"), method=c("permutation"), nsim = 1000)
cospeciation(lyve_tree,enterobase_tree, distance = c("RF"), method=c("permutation"), nsim = 1000)
cospeciation(cfsan_tree, enterobase_tree, distance = c("RF"), method=c("permutation"), nsim = 1000)

# Test using Subtree pruning and regrafting (SPR)
cospeciation(ksnp_tree, cfsan_tree, distance = c("SPR"), method=c("permutation"), nsim = 1000)
cospeciation(ksnp_tree, enterobase_tree, distance = c("SPR"), method=c("permutation"), nsim = 1000)
cospeciation(lyve_tree,ksnp_tree, distance = c("SPR"), method=c("permutation"), nsim = 1000)
cospeciation(lyve_tree,cfsan_tree, distance = c("SPR"), method=c("permutation"), nsim = 1000)
cospeciation(lyve_tree,enterobase_tree, distance = c("SPR"), method=c("permutation"), nsim = 1000)
cospeciation(cfsan_tree, enterobase_tree, distance = c("SPR"), method=c("permutation"), nsim = 1000)

# Example plot of cospeciation results
plot(cospeciation(ksnp_tree, enterobase_tree, distance = c("RF")))


#
## Compare trees with all.equal.phylo
#
all.equal.phylo(lyve_tree,ksnp_tree)
all.equal.phylo(lyve_tree,cfsan_tree)
all.equal.phylo(lyve_tree,enterobase_tree)
all.equal.phylo(ksnp_tree, cfsan_tree)
all.equal.phylo(ksnp_tree, enterobase_tree)
all.equal.phylo(cfsan_tree, enterobase_tree)

#
## Plots - 
#
comparePhylo(lyve_tree,ksnp_tree, plot=TRUE)
comparePhylo(lyve_tree,cfsan_tree, plot=TRUE)
comparePhylo(lyve_tree,enterobase_tree, plot=TRUE)
comparePhylo(ksnp_tree, cfsan_tree, plot=TRUE)
comparePhylo(ksnp_tree, enterobase_tree, plot=TRUE)
comparePhylo(cfsan_tree, enterobase_tree, plot=TRUE)

# 
## Plots -line connecting nodes between trees
# http://phytools.org/mexico2018/ex/12/Plotting-methods.html
# Compare trees with all.equal.phylo
plot(cophylo(lyve_tree,ksnp_tree))
plot(cophylo(lyve_tree,cfsan_tree))
plot(cophylo(lyve_tree,enterobase_tree))
plot(cophylo(ksnp_tree, cfsan_tree))
plot(cophylo(ksnp_tree, enterobase_tree))
plot(cophylo(cfsan_tree, enterobase_tree))

# Side-by-side comparison with phylogenetic trees
obj <- cophylo(ksnp_tree, cfsan_tree, print= TRUE)
plot(obj,link.type="curved",link.lwd=3,link.lty="solid",
     link.col="grey",fsize=0.8)
nodelabels.cophylo(which="left",frame="circle",cex=0.8)
nodelabels.cophylo(which="right",frame="circle",cex=0.8)

#
##
###
#### Combine rooted and cleaned trees
###
##
#
combined_rooted_trees <- c(ksnp_tree_rooted,enterobase_tree_rooted,cfsan_tree_rooted, lyve_tree_rooted)
combined_cleaned_trees <- c(ksnp_tree,enterobase_tree,cfsan_tree, lyve_tree)
names(combined_cleaned_trees) <- c("ksnp","enterobase","cfsan","lyveset")

densityTree(combined_rooted_trees,type="cladogram",nodes="intermediate")
densityTree(combined_cleaned_trees,type="cladogram",nodes="intermediate")
densityTree(combined_rooted_trees,use.edge.length=FALSE,type="phylogram",nodes="inner", alpha = 0.3)

# Load updated metadata file 
#SRA_metadata <- read.csv("SRA_present.csv", header = FALSE, stringsAsFactors = FALSE)
# Calculate related tree distance
#relatedTreeDist(combined_trees, as.data.frame(SRA_metadata), checkTrees = TRUE)

#write.csv(lyve_tree$tip.label, "lyve_tree_nodes.csv")
#write.csv(ksnp_tree$tip.label, "ksnp3_nodes.csv")

# png(filename = "List_NY_lyveset_tree.png", res = 300,width = 800, height = 800)
# plotTree(lyve_tree, label.offset =1)





###
####
#######  Treespace
####
###
# https://cran.r-project.org/web/packages/treespace/vignettes/introduction.html

combined_treespace <- treespace(combined_rooted_trees, nf=3) # , return.tree.vectors = TRUE

#test <- as.treeshape(dataset1_tree_vector)
table.image(combined_treespace$D)
table.value(combined_treespace$D, nclass=5, method="color", symbol="circle", col=redpal(6))
plotGroves(combined_treespace$pco,lab.show=TRUE, lab.cex=1.5)
combined_treespace_groves <- findGroves(combined_treespace)
plotGrovesD3(combined_treespace_groves)


#aldous.test(combined_rooted_trees)

colless.test(combined_treespace_groves, alternative="greater")
likelihood.test(combined_treespace, alternative='greater')


##
###
#### ggtree
###
##

# Mutate to create new column with selected outbreak group
SRA_metadata <- as_tibble(SRA_metadata)
SRA_metadata <- SRA_metadata %>% mutate(Group = ifelse(SNP.cluster == "PDS000000366.382" , "Outbreak", "Other"))








