library(ggplot2)
library(reshape2)

IGHV <- read.table("result/PotentialHV.tsv")
IGHV <- data.frame(start = IGHV$sstart, end = IGHV$send, group = "IGHV")

IGHD <- data.frame(t(data.frame(list(c(2083,2151),
          c( 3162,3245),
          c( 4132,4214),
          c( 5282,5367),
          c( 6190,6282)))))
colnames(IGHD) <- c("start", "end")
IGHD <- IGHD + 534854
IGHD$group = "IGHD"
IGHJ <- data.frame(start = c(541833), end = c(541928), group = "IGHJ")

IGH <- rbind(IGHV, IGHD, IGHJ)
IGH[c('start', 'end')] <- lapply(IGH[c('start', 'end')], as.numeric)

class(IGH$start)

ggplot(IGH, aes(x = start, xend = end, y = 1, yend = 1, color = group)) +
  geom_segment(size = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(0, 541928, by = 5000)) +
  labs(title = "IGH gene segments",
       x = "Genomic position",
       y = "Gene segment") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("plot/IGH_gene_segments.png", width = 10, height = 4)
