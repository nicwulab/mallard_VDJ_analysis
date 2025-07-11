library(ggplot2)
library(stringr)
library(patchwork)  

PlotSum <- function(fr_tb, fr_id){
  P1 <- ggplot(fr_tb, aes(x=sacc, fill = sample)) + geom_bar() + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    ggtitle("Counts") + xlab("Contig ID") + ylab("Count") + theme(legend.position = "left")
  P2 <- ggplot(fr_tb, aes(x=sacc, y = pident, fill = sample)) + geom_violin() + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    ggtitle("Identity") + xlab("Contig ID") + ylab("Percent Identity") + theme(legend.position = "none") 
  P3 <- ggplot(fr_tb, aes(x=sacc, y = length, fill = sample)) + geom_violin() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    ggtitle("Length") + xlab("Contig ID") + ylab("Length") + theme(legend.position = "none")
  P4 <- ggplot(fr_tb, aes(x=sacc, y = mismatch, fill = sample)) + geom_violin() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    ggtitle("Mismatch") + xlab("Contig ID") + ylab("mismatch") + theme(legend.position = "none")
  P5 <- ggplot(fr_tb, aes(x=sacc, y = gaps, fill = sample)) + geom_violin() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    ggtitle("Gaps") + xlab("Contig ID") + ylab("gaps") + theme(legend.position = "none")
  GGlayout <- 'ABBCC
  ADDEE
  '
  P1 + P2 + P3 + P4 + P5 + plot_layout(design = GGlayout)
  ggsave(paste("plot/Blast_summary_", fr_id, ".png", sep=""), width=7.5, height=5.3)
}


head(fr_tb)


tb_match_lco <- data.frame()

fr_files <- paste("result/", grep("fr..tsv", list.files('result'), value=TRUE), sep = "")
for(fr_loc in fr_files){
  fr_id <- substr(fr_loc, 8, 10)
  fr_tb <- read.table(fr_loc, header=TRUE, sep="\t") 
  colnames(fr_tb) <- c('qacc', 'sacc', 'pident', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'length', 'mismatch', 'gaps')
  fr_tb$sample <- as.data.frame(str_split_fixed(fr_tb$qacc, ":", 2))[[1]]
  fr_tb$seg <- fr_id
  tb_match_lco <- rbind(tb_match_lco, fr_tb)
  PlotSum(fr_tb, fr_id)
}

fr_files <- paste("result/", grep("cdr..tsv", list.files('result'), value=TRUE), sep = "")
for(fr_loc in fr_files){
  fr_id <- substr(fr_loc, 8, 11)
  fr_tb <- read.table(fr_loc, header=TRUE, sep="\t")
  colnames(fr_tb) <- c('qacc', 'sacc', 'pident', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'length', 'mismatch', 'gaps')
  fr_tb$sample <- as.data.frame(str_split_fixed(fr_tb$qacc, ":", 2))[[1]]
  fr_tb$seg <- fr_id
  tb_match_lco <- rbind(tb_match_lco, fr_tb)
  PlotSum(fr_tb, fr_id)
}

head(tb_match_lco)

P1 <- ggplot(tb_match_lco, aes(x = sstart, y = seg, color = seg)) + geom_point() + 
  facet_wrap(~sacc, scales = "free") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

ggsave("plot/SegAlignment2.png", width=10, height=4.5)


##

tb_v <- read.table("result/potential_fr_complete.tsv", sep = "\t", header = TRUE)[1:10] 
tb_v$sfd <- as.data.frame(str_split_fixed(tb_v$X, "_", 2))[[1]]


ggplot(tb_v[tb_v$sfd == 'ptg000189l',], aes(fr1_s, sfd))+
  geom_segment(aes(x = fr1_s, xend = fr3_e, y = sfd, yend = sfd, size =5), color = 'black') +
  geom_segment(aes(x = fr4_s, xend = fr4_e, y = sfd, yend = sfd, size =5), color = 'salmon') +
  theme_bw()

head(tb_v)



# plot the mapped contigs


