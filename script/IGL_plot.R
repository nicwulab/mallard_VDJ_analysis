library(ggplot2)
library(reshape2)

IGLV <- read.table("result/PotentialLV.tsv")
IGLV <- data.frame(start = IGLV$sstart, end = IGLV$send, group = "IGLV", has_sig = IGLV$has_sig)

IGLJ <- data.frame(start = c(9508724), end = c(9508789), group = "IGLJ", has_sig = "True")

IGL <- rbind(IGLV, IGLJ)
IGL[c('start', 'end')] <- lapply(IGL[c('start', 'end')], as.numeric)

class(IGL$start)

ggplot(IGL, aes(x = start, xend = end, y = 1, yend = 1, color = group)) +
  geom_segment(size = 5) +
  geom_text(data = IGL[IGL$has_sig == "True",], aes(x = (start + end)/2, y = 1, label = "*"), size = 5, vjust = .8, color = 'black') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(0, max(IGL$end) , by = 5000)) +
  labs(title = "IGL gene segments",
       x = "Genomic position",
       y = "Gene segment") +
  theme(plot.title = element_text(hjust = 0.5))

