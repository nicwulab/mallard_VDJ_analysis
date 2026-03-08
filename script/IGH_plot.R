library(ggplot2)
library(stringr)
library(reshape2)


IGHV <- grep(">", read.table('final/IGHV.fa', header = FALSE, sep= '\t')[[1]], value = TRUE)
IGHD <- grep(">", read.table('final/IGHD.fa', header = FALSE, sep= '\t')[[1]], value = TRUE)
IGHJ <- grep(">", read.table('final/IGHJ.fa', header = FALSE, sep= '\t')[[1]], value = TRUE)

IGH <- c(IGHV, IGHD, IGHJ)
IGH.tb <- data.frame(str_split_fixed(IGH, ' ', 5))
colnames(IGH.tb) <- c('ID', 'Chr', 'sstart', 'send', 'Direction')
IGH.tb$Direction[which(IGH.tb$Direction == "")] <- "+"

IGH.tb$ID <- str_replace(IGH.tb$ID, ">", "")
# get the first 4 letters 
IGH.tb$Type <- str_sub(IGH.tb$ID, 1, 4)

IGH.tb$sstart <- as.numeric(IGH.tb$sstart)
IGH.tb$send <- as.numeric(IGH.tb$send)

IGH.tb$Fct <- 'Pseudogene'
IGH.tb$Fct[grep("IGHV", IGH.tb$ID, invert = TRUE)] <- 'Functional'
IGH.tb$Fct[1] <- 'Functional'


ggplot(IGH.tb, aes(x = sstart, xend = send, y = 1, yend = 1, color = Type)) +
  geom_line(color = 'black', size = 0.5) +
  geom_segment(fill ='black', size = 20) +
  geom_segment(data = IGH.tb[IGH.tb$Direction == "-", ],
    aes(x = send, xend = sstart, y = 1.03, yend = 1.03), color = 'black',
    size = 1, arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(aes(x = (sstart + send)/2, y = 1.04, label = ID), angle = 90,
  vjust = 0.5, fontface = "italic", hjust=0, size = 3) +
  geom_label(data = IGH.tb[IGH.tb$Fct == 'Functional',], fontface = "italic",
    aes(x = (sstart + send)/2, y = 1.04, label = ID), angle = 90, vjust = 0.5, hjust=0, size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(min(IGH.tb$sstart), max(IGH.tb$sstart), by = 5000)) +
  labs(title = "IGH gene segments",
       x = "Genomic position",
       y = "Gene segment") +
  theme(plot.title = element_text(hjust = 0.5),
      legend.position = 'none') +
  coord_cartesian(ylim = c(.98, 1.1)) +
  scale_color_manual(values = c("IGHV" = "#0B2D72", "IGHD" = "#40513B", "IGHJ" = "#740A03"))

ggsave("plot/IGH_gene_segments_no_direction.png", width = 19.4, height = 2.68)


# IGL
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")

# Load the library
library(Biostrings)

fasta_file <- "final/IGLV.pro"
sequences <- readAAStringSet(fasta_file)

sequences_with_stop <- grep("\\*", as.character(sequences))
WithStop <- str_split_fixed(names(sequences), ' ', 2)[,1][sequences_with_stop]

IGLV <- grep(">", read.table('final/IGLV.fa', header = FALSE, sep= '\t')[[1]], value = TRUE)
IGLJ <- grep(">", read.table('final/IGLJ.fa', header = FALSE, sep= '\t')[[1]], value = TRUE)
IGL <- c(IGLV, IGLJ)

IGL.tb <- data.frame(str_split_fixed(IGL, ' ', 8))
colnames(IGL.tb) <- c('ID', 'Chr', 'sstart', 'send', 'Direction', 'RSS', 'Spacer', 'Distance')
IGL.tb$Direction[which(IGL.tb$Direction == "")] <- "+"
IGL.tb$ID <- str_replace(IGL.tb$ID, ">", "")
# get the first 4 letters
IGL.tb$Type <- str_sub(IGL.tb$ID, 1, 4)
IGL.tb$sstart <- as.numeric(IGL.tb$sstart)
IGL.tb$send <- as.numeric(IGL.tb$send)
IGL.tb$Fct <- 'Pseudogene'
IGL.tb$Distance <- as.numeric(IGL.tb$Distance)

IGL.tb$Fct[IGL.tb$Distance<=20] <- 'Functional'
IGL.tb$Fct[IGL.tb$Type=='IGLJ'] <- 'Functional'

IGL.tb$Fct[which(IGL.tb$ID %in% WithStop)] <- 'Pseudogene'

# get last 9 characters of RSS
IGL.tb$RSS_pattern <- paste0(str_sub(IGL.tb$RSS,1,7), "-",IGL.tb$Spacer ,"-", str_sub(IGL.tb$RSS, -9))
IGL.tb$RSS_pattern[IGL.tb$RSS_pattern=='--'] = ''

IGL.up <- IGL.tb[seq(1, nrow(IGL.tb), 2)-1,]
IGL.down <- IGL.tb[seq(1, nrow(IGL.tb), 2),]

flipList <- c("IGLV1-70", "IGLV1-18")
IGL.down <- rbind(IGL.down, IGL.up[IGL.up$ID %in% flipList,])
IGL.up <- IGL.up[!IGL.up$ID %in% flipList,]

# remove functional gene from IG.up
IGL.down <- IGL.down[!IGL.down$Fct == 'Functional',]

ggplot(IGL.tb, aes(x = sstart, xend = send, y = 1, yend = 1, color = Type)) +
  geom_line(color = 'black', size = 0.5) +
  geom_segment(fill ='black', size = 20) +
  geom_segment(data = IGL.tb[IGL.tb$Direction == "-", ],
    aes(x = send, xend = sstart, y = 1.03, yend = 1.03), color = 'black',
    size = 1, arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = IGL.up,
    aes(x = (sstart + send)/2, y = 1.04, label = ID), angle = 90, 
    fontface = "italic", vjust = 0.5, hjust=0, size = 3) +
  geom_text(data = IGL.down,
    aes(x = (sstart + send)/2, y = .96, label = ID), angle = 90,
    fontface = "italic", vjust = 0.5, hjust=1, size = 3) +
  geom_label(data = IGL.tb[IGL.tb$Fct == 'Functional',], fontface = "italic",
    aes(x = (sstart + send)/2, y = 1.04, label = ID), angle = 90, vjust = 0.5, hjust=0, size = 3) +
  geom_text(aes(x = (sstart + send)/2,
    y = 1.11, label = RSS_pattern), angle = 90, vjust = 0.5, hjust=0, size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(min(IGL.tb$sstart), max(IGL.tb$send), by = 5000)) +
  labs(title = "IGL gene segments",
       x = "Genomic position",
       y = "Gene segment") +
  theme(plot.title = element_text(hjust = 0.5),
      legend.position = 'none') +
  coord_cartesian(ylim = c(.90, 1.3)) +
  scale_color_manual(values = c("IGLV" = "#0B2D72", "IGLJ" = "#740A03"))

ggsave("plot/IGL_gene_segments_with_RSS.png", width = 19, height = 5.13)
