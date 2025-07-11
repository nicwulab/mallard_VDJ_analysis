library(ggtree)

tree <- read.tree("result/Rag1_Uniprot.dnd")

ggplot(tree, aes(x, y)) + geom_tree() + 
  theme_tree() + geom_tiplab(size=5, color="purple") +
  xlim(NA, 1)
