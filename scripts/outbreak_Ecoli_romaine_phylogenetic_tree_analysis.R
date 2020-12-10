
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
library(ape)
library(ggplot2)
library(tidyverse)

# Set seed for reproducibility 
set.seed(23)

# Load metadata file
SRA_metadata <- read.table("Pipeline_results/Ecoli_romaine_outbreak/exported_trees/Ecoli_romaine_metadata.csv",header = TRUE, fill = FALSE ,sep = ',',stringsAsFactors = TRUE)
# typeof(SRA_metadata)
# View(SRA_metadata)
# typeof(as.data.frame(SRA_metadata))

# Read in phylogenetic trees
lyve_tree <- read.tree(file = "Pipeline_results/Ecoli_romaine_outbreak/exported_trees/lyveset_NJ.newick")
# kSNP3 tree.NJ.tre, tree.ML.tre, tree.core.tre, tree.parsimony.tre
ksnp_tree <- read.tree(file = "Pipeline_results/Ecoli_romaine_outbreak/exported_trees/ksnp3_NJ.newick")
# Cfsan
cfsan_tree <- read.tree(file = "Pipeline_results/Ecoli_romaine_outbreak/exported_trees/cfsan_NJ.newick")
# Enterobase
enterobase_tree <- read.tree(file = "Pipeline_results/Ecoli_romaine_outbreak/exported_trees/enterobase_NJ.newick")



# Combine trees
combined_trees <- c(lyve_tree,ksnp_tree,cfsan_tree,enterobase_tree)
# Combine trees from single dataset into vector
dataset1_tree_vector <- c(lyve_tree,ksnp_tree,cfsan_tree,enterobase_tree)
dataset1_tree_vector <- c(as.phylo(lyve_tree),as.phylo(ksnp_tree),as.phylo(cfsan_tree),as.phylo(enterobase_tree))

#
##
### Code for subsetting trees with unmatched nodes
## Need to automate this 
# SetDiff
setdiff(cfsan_tree$tip.label, ksnp_tree$tip.label)

## Check for sample matches
# Find samples not in cfsan_snp (lowest number of tips)
all_SRA_to_drop = c()

# SRA_to_drop <- unique(enterobase_tree$tip.label[! enterobase_tree$tip.label %in% cfsan_tree$tip.label])
# all_SRA_to_drop = c(all_SRA_to_drop,SRA_to_drop)
SRA_to_drop <- unique(enterobase_tree$tip.label[! enterobase_tree$tip.label %in% ksnp_tree$tip.label])
all_SRA_to_drop = c(all_SRA_to_drop,SRA_to_drop)
SRA_to_drop <- unique(cfsan_tree$tip.label[! cfsan_tree$tip.label %in% enterobase_tree$tip.label])
all_SRA_to_drop = c(all_SRA_to_drop,SRA_to_drop)

# SRA_to_drop <- unique(lyve_tree$tip.label[! lyve_tree$tip.label %in% ksnp_tree$tip.label])
# SRA_to_drop <- unique(ksnp_tree$tip.label[! ksnp_tree$tip.label %in% lyve_tree$tip.label])

all_SRA_to_drop <- unique(all_SRA_to_drop)
#SRA_to_drop <- unique(cfsan_tree$tip.label[! cfsan_tree$tip.label %in% ksnp_tree$tip.label])
#SRA_to_drop <- unique(ksnp_tree$tip.label[! ksnp_tree$tip.label %in% cfsan_tree$tip.label])
#SRA_to_drop <- unique(cfsan_tree$tip.label[! cfsan_tree$tip.label %in% lyve_tree$tip.label])


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

# png(filename = "Ecoli_romaine_lyveset_tree.png", res = 300,width = 800, height = 800)
# plotTree(lyve_tree, label.offset =1)







##
###
#### ggtree
###
##

# Mutate to create new column with selected outbreak group
SRA_metadata <- as_tibble(SRA_metadata)
SRA_metadata <- SRA_metadata %>% mutate(Group = ifelse(SNP.cluster == "PDS000046273.15" , "Outbreak", "Other"))
SRA_metadata$Group <- SRA_metadata$Newick_label

# Lyveset
#ggtree(lyve_tree,branch.length='none')  + theme_tree2() + geom_nodepoint(color="#b5e521", alpha=1/4, size=10) + geom_nodelab(geom = "text")


# lyve
plyve <- ggtree(lyve_tree) %<+% SRA_metadata 
plyve2 <- plyve + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint( size = 2) +  
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))
plyve2
ggsave("Pipeline_results/Ecoli_romaine_outbreak/lyve_tree_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()

# lyve no lengths
plyve <- ggtree(lyve_tree, branch.length = "none") %<+% SRA_metadata 
plyve2 <- plyve + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint( size = 2) + 
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))
plyve2
ggsave("Pipeline_results/Ecoli_romaine_outbreak/lyve_tree_nolengths_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()


plyve2 <- plyve + geom_tiplab(offset = .6, hjust = .5) + geom_text(aes(x=branch, label = branch))
plyve2

ggtree(ksnp_tree) + theme_tree2()

# cfsan
pcfsan <- ggtree(cfsan_tree) %<+% SRA_metadata 
pcfsan2 <- pcfsan + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint( size = 2) + 
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))
pcfsan2
ggsave("Pipeline_results/Ecoli_romaine_outbreak/cfsan_tree_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()

# cfsan_nolengths
pcfsan <- ggtree(cfsan_tree, branch.length = "none") %<+% SRA_metadata 
pcfsan2 <- pcfsan + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint( size = 2) +  
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))
pcfsan2
ggsave("Pipeline_results/Ecoli_romaine_outbreak/cfsan_tree_nolengths_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()


# enterobase
penterobase <- ggtree(enterobase_tree) %<+% SRA_metadata 
penterobase2 <- penterobase + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint( size = 2) + 
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))
penterobase2
ggsave("Pipeline_results/Ecoli_romaine_outbreak/enterobase_tree_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()

# enterobase_nolengths
penterobase <- ggtree(enterobase_tree, branch.length = "none") %<+% SRA_metadata 
penterobase2 <- penterobase + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint( size = 2) + 
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))
penterobase2
ggsave("Pipeline_results/Ecoli_romaine_outbreak/enterobase_tree_nolengths_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()


# ksnp
pksnp <- ggtree(ksnp_tree) %<+% SRA_metadata 
pksnp2 <- pksnp + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint( size = 2) +  
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))
pksnp2
ggsave("Pipeline_results/Ecoli_romaine_outbreak/ksnp_tree_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()

# ksnp_nolengths
pksnp <- ggtree(ksnp_tree, branch.length = "none") %<+% SRA_metadata 
pksnp2 <- pksnp + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint( size = 2) + 
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))
pksnp2
ggsave("Pipeline_results/Ecoli_romaine_outbreak/ksnp_tree_nolengths_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()


##
###
#### ggtree combined tree figures
###
##


pfacet_tree <- ggtree(combined_cleaned_trees) + facet_wrap( ~.id, scale="free") + theme_tree2()
pfacet_tree
ggsave("Pipeline_results/Ecoli_romaine_outbreak/combined_trees_wblengths_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()

# no branch lengths
pfacet_tree2 <- ggtree(combined_cleaned_trees,branch.length='none') + facet_wrap( ~.id, scale="free") + theme_tree2()
pfacet_tree2
ggsave("Pipeline_results/Ecoli_romaine_outbreak/combined_trees_nolengths_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()


pdense_tree <- ggdensitree(combined_cleaned_trees, alpha=.3, colour='steelblue') + 
  geom_tiplab(size=3) + xlim(0, 45)
pdense_tree
ggsave("Pipeline_results/Ecoli_romaine_outbreak/combined_trees_dense_tree_Ecoli_romaine.png", width = 50, height = 80, units = "cm")
dev.off()



#
##
###
#### Info on tree
###
##
#
summary(lyve_tree)
sum(lyve_tree$edge.length)

summary(cfsan_tree)
sum(cfsan_tree$edge.length)

summary(enterobase_tree)
sum(enterobase_tree$edge.length)

summary(ksnp_tree)
sum(ksnp_tree$edge.length)


# Calculate co-speciation (RF distance) between trees
# Make function to automate these pairwise comparisons within a vector of trees
# Save results in table
# Compare trees with cospeciation
cospeciation(ksnp_tree, cfsan_tree, distance = c("RF","SPR"))
cospeciation(ksnp_tree, enterobase_tree, distance = c("RF","SPR"))
cospeciation(lyve_tree,ksnp_tree, distance = c("RF","SPR"))
cospeciation(lyve_tree,cfsan_tree, distance = c("RF","SPR"))
cospeciation(lyve_tree,enterobase_tree, distance = c("RF","SPR"))
plot(cospeciation(ksnp_tree, enterobase_tree, distance = c("RF")))
cospeciation(cfsan_tree, enterobase_tree, distance = c("RF","SPR"))


# Compare trees with all.equal.phylo
all.equal.phylo(lyve_tree,ksnp_tree)
all.equal.phylo(lyve_tree,cfsan_tree)
all.equal.phylo(lyve_tree,enterobase_tree)
all.equal.phylo(ksnp_tree, cfsan_tree)
all.equal.phylo(ksnp_tree, enterobase_tree)
all.equal.phylo(cfsan_tree, enterobase_tree)


# Plots 
comparePhylo(lyve_tree,ksnp_tree, plot=TRUE)
comparePhylo(lyve_tree,cfsan_tree, plot=TRUE)
comparePhylo(lyve_tree,enterobase_tree, plot=TRUE)
comparePhylo(ksnp_tree, cfsan_tree, plot=TRUE)
comparePhylo(ksnp_tree, enterobase_tree, plot=TRUE)
comparePhylo(cfsan_tree, enterobase_tree, plot=TRUE)

# Side-by-side tree node comparison plots
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




#
##
### Test OTU grouping
##
#
outbreak_group <- SRA_metadata %>%
  filter(SNP.cluster == "PDS000000366.382")  %>%
  select(Newick_label)

#Just the sample Ids
lyve_tree_w_meta <- groupOTU(lyve_tree ,outbreak_group, group_name = "Outbreak")

p <- ggtree(lyve_tree_w_meta, aes(color=Outbreak)) +
  scale_color_manual(values = c("#efad29", "#63bbd4")) +
  geom_nodepoint(color="black", size=0.1) +
  geom_tiplab(size=2, color="black")
p



##
###
#### Extra code
###
##

# Branch lengths
plot(compute.brlen(lyve_tree))
# Branch times
#This function computes the branch lengths of a tree giving its branching times (aka node ages or heights).
plot(compute.brtime(lyve_tree))

# dist.topo
dist.topo(combined_cleaned_trees)

dnd1 <- as.dendrogram(lyve_tree)
dnd2 <- as.dendrogram(ksnp_tree)

plotTree(enterobase_tree)
