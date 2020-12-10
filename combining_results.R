lyveset_both_trees <- c(combined_cleaned_trees[[4]], core_combined_cleaned_trees[[4]])
names(lyveset_both_trees) <- c("default lyveset", "core lyveset")


pfacet_tree <- ggtree(lyveset_both_trees,) + facet_wrap( ~.id, scale="free") + theme_tree2()
pfacet_tree
ggsave("Pipeline_results/Listeria_NY_and_outbreak_geography/lyveset_combined_trees_wblengths_Listeria_NY_and_outbreak.png", width = 50, height = 80, units = "cm")
dev.off()



pfacet_tree <- ggtree(lyveset_both_trees,branch.length = "none") + facet_wrap( ~.id, scale="free") + theme_tree2()
pfacet_tree
ggsave("Pipeline_results/Listeria_NY_and_outbreak_geography/lyveset_combined_trees_nolengths_Listeria_NY_and_outbreak.png", width = 50, height = 80, units = "cm")
dev.off()
