url <- paste0("https://raw.githubusercontent.com/TreeViz/",
            "metastyle/master/design/viz_targets_exercise/")
x <- read.tree(paste0(url, "tree_boots.nwk"))
info <- read.csv(paste0(url, "tip_data.csv"))

p <- ggtree(x) %<+% info + xlim(-.1, 4)
p2 <- p + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint(aes(shape = trophic_habit, color = trophic_habit, size = mass_in_kg)) + 
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))



d2 <- read.csv(paste0(url, "inode_data.csv"))

p2 %<+% d2 + geom_label(aes(label = vernacularName.y, fill = posterior)) + 
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "YlGnBu"))
