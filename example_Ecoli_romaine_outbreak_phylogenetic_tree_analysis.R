# Required R packages
library(ape)
library(phytools)
install.packages('TreeDist')
library(TreeDist)
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

# Add root to tree
lyve_tree_rooted <- root(lyve_tree,1, r = TRUE)
ksnp_tree_rooted <- root(ksnp_tree,1, r = TRUE)
cfsan_tree_rooted <- root(cfsan_tree,1, r = TRUE)
enterobase_tree_rooted <- root(enterobase_tree,1, r = TRUE)


combined_trees_clean <- c(lyve_tree,ksnp_tree,cfsan_tree,enterobase_tree)



#
##
### TreeDist
## Generalized Robinson-Foulds distance
#
VisualizeMatching(SharedPhylogeneticInfo, lyve_tree, ksnp_tree, 
                  Plot = TreeDistPlot, matchZeros = FALSE)


SharedPhylogeneticInfo(lyve_tree, ksnp_tree)
MutualClusteringInfo(lyve_tree, ksnp_tree)
NyeSimilarity(lyve_tree, ksnp_tree)
JaccardRobinsonFoulds(lyve_tree, ksnp_tree)
MatchingSplitDistance(lyve_tree, ksnp_tree)
MatchingSplitInfoDistance(lyve_tree, ksnp_tree)


VisualizeMatching(JaccardRobinsonFoulds, lyve_tree, ksnp_tree, 
                  Plot = TreeDistPlot, matchZeros = FALSE)

#
##
### TreeDist
## Using a suitable distance metric, projecting distances
#

# Tree colors
library('TreeTools', quietly = TRUE, warn.conflicts = FALSE)
treeNumbers <- c(1:4)
spectrum <- viridisLite::plasma(4)
treeCols <- spectrum[treeNumbers]

# calculate distances
distances <- ClusteringInfoDistance(combined_trees_clean)
distances <- RobinsonFoulds(combined_trees_clean)


distances <- as.dist(Quartet::QuartetDivergence(Quartet::ManyToManyQuartetAgreement(combined_trees_clean), similarity = FALSE))



# Projecting distances
#Then we need to reduce the dimensionality of these distances. We’ll start out with a 12-dimensional projection; if needed, we can always drop higher dimensions.

#Principal components analysis is quick and performs very well:
  
projection <- cmdscale(distances, k = 3)
# Alternative projection methods do exist, and sometimes give slightly better projections. isoMDS() performs non-metric multidimensional scaling (MDS) with the Kruskal-1 stress function (Kruskal, 1964):
kruskal <- MASS::isoMDS(distances, k = 3)
projection <- kruskal$points
#whereas sammon(), one of many metric MDS methods, uses Sammon’s stress function (Sammon, 1969):
  
sammon <- MASS::sammon(distances, k = 3)
projection <- sammon$points
#That’s a good start. It is tempting to plot the first two dimensions arising from this projection and be done:
  
par(mar = rep(0, 4))
plot(projection,
     asp = 1, # Preserve aspect ratio - do not distort distances
     ann = FALSE, axes = FALSE, # Don't label axes: dimensions are meaningless
     col = treeCols, pch = 16
)

#
## Identifying clusters
#


# A quick visual inspection suggests at least two clusters, with the possibility of further subdivision
# of the brighter trees. But visual inspection can be highly misleading (Smith, 2021). 
# We must take a statistical approach. A combination of partitioning around medoids and hierarchical 
# clustering with minimax linkage will typically find a clustering solution that is close to optimal, 
# if one exists (Smith, 2021).
library(protoclust)

possibleClusters <- 3:10

# Had to choose static K value
pamClusters <- lapply(possibleClusters, function (x) cluster::pam(distances, k = 3))

pamSils <- vapply(pamClusters, function (pamCluster) {
  mean(cluster::silhouette(pamCluster)[, 3])
}, double(1))

bestPam <- which.max(pamSils)
pamSil <- pamSils[bestPam]
pamCluster <- pamClusters[[bestPam]]$cluster

hTree <- protoclust::protoclust(distances)
hClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = 3))
hSils <- vapply(hClusters, function (hCluster) {
  mean(cluster::silhouette(hCluster, distances)[, 3])
}, double(1))


bestH <- which.max(hSils)
hSil <- hSils[bestH]
hCluster <- hClusters[[bestH]]

plot(pamSils ~ possibleClusters,
     xlab = 'Number of clusters', ylab = 'Silhouette coefficient',
     ylim = range(c(pamSils, hSils)))
points(hSils ~ possibleClusters, pch = 2)
legend('topright', c('PAM', 'Hierarchical'), pch = 1:2)

# Silhouette coefficients of < 0.25 suggest that structure is not meaningful; > 0.5 denotes good evidence 
# of clustering, and > 0.7 strong evidence (Kaufman & Rousseeuw, 1990). The evidence for the visually 
# apparent clustering is not as strong as it first appears. Let’s explore our two-cluster hierarchical
# clustering solution anyway.

cluster <- hClusters[[2 - 1]]
#We can visualize the clustering solution as a tree:
  
class(hTree) <- 'hclust'
par(mar = c(0, 0, 0, 0))
plot(hTree, labels = FALSE, main = '')
points(seq_along(trees), rep(1, length(trees)), pch = 16,
       col = spectrum[hTree$order])


#Another thing we may wish to do is to take the consensus of each cluster:
  
par(mfrow = c(1, 2), mar = rep(0.2, 4))
col1 <- spectrum[mean(treeNumbers[cluster == 1])]
col2 <- spectrum[mean(treeNumbers[cluster == 2])]
plot(consensus(trees[cluster == 1]), edge.color = col1, edge.width = 2, tip.color = col1)
plot(consensus(trees[cluster == 2]), edge.color = col2, edge.width = 2, tip.color = col2)



# Validating a projection
# Now let’s evaluate whether our plot of tree space is representative. First we want to know how many dimensions are necessary to adequately represent the true distances between trees. We hope for a trustworthiness × continuity score of > 0.9 for a usable projection, or > 0.95 for a good one.
library(TreeTools)
# ProjectionQuality doesn't work with regular TreeDist
#remotes::install_github("ms609/TreeDist")

txc <- vapply(1:12, function (k) {
  newDist <- dist(projection[, seq_len(k)])
  TreeTools::ProjectionQuality(distances, newDist, 10)['TxC']
}, 0)
plot(txc, xlab = 'Dimension')
abline(h = 0.9, lty = 2)


# To help establish visually what structures are more likely to be genuine, we might also choose to calculate a minimum spanning tree:
mstEnds <- MSTEdges(distances)

# Let’s plot the first five dimensions of our tree space, highlighting the convex hulls of our clusters:
plotSeq <- matrix(0, 5, 5)
plotSeq[upper.tri(plotSeq)] <- seq_len(5 * (5 - 1) / 2)
plotSeq <- t(plotSeq[-5, -1])
plotSeq[c(5, 10, 15)] <- 11:13
layout(plotSeq)
par(mar = rep(0.1, 4))

for (i in 2:4) for (j in seq_len(i - 1)) {
  # Set up blank plot
  plot(projection[, j], projection[, i], ann = FALSE, axes = FALSE, frame.plot = TRUE,
       type = 'n', asp = 1, xlim = range(projection), ylim = range(projection))
  
  # Plot MST
  apply(mstEnds, 1, function (segment)
    lines(projection[segment, j], projection[segment, i], col = "#bbbbbb", lty = 1))
  
  # Add points
  points(projection[, j], projection[, i], pch = 16, col = treeCols)
  
  # Mark clusters
  for (clI in unique(cluster)) {
    inCluster <- cluster == clI
    clusterX <- projection[inCluster, j]
    clusterY <- projection[inCluster, i]
    hull <- chull(clusterX, clusterY)
    polygon(clusterX[hull], clusterY[hull], lty = 1, lwd = 2,
            border = '#54de25bb')
  }
}
# Annotate dimensions
plot(0, 0, type = 'n', ann = FALSE, axes = FALSE)
text(0, 0, 'Dimension 2')
plot(0, 0, type = 'n', ann = FALSE, axes = FALSE)
text(0, 0, 'Dimension 3')
plot(0, 0, type = 'n', ann = FALSE, axes = FALSE)
text(0, 0, 'Dimension 4')



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
plot(cospeciation(ksnp_tree, enterobase_tree, distance = c("RF"),method=c("permutation"), nsim = 1000))


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
library(treespace)


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
### Test OTU grouping - not fully working yet.
##
# Enrique working on this

# Mutate to create new column with selected outbreak group
SRA_metadata <- as_tibble(SRA_metadata)
SRA_metadata <- SRA_metadata %>% mutate(Group = ifelse(SNP.cluster == "PDS000000366.382" , "Outbreak", "Other"))


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







