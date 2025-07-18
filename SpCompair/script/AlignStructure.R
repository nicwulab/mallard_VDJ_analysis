library(ggplot2)
library(reshape2)

TB <- read.csv('result/Aligned_SecStruct.csv')

Odds <- c()
NonOdds <- c()
for(i in c(2:ncol(TB))){
  if(i%%2 == 0){
    Odds <- c(Odds, i)
  } else {
    NonOdds <- c(NonOdds, i)
  }
}

TB$human <- row.names(TB)

TB.resi <- melt(TB[c(1, Odds)], id.vars = 'human')
TB.Stru <- melt(TB[c(1, NonOdds)], id.vars = 'human')

TB.Stru$variable <- factor(TB.Stru$variable, 
                          levels = sort(as.character(unique(TB.Stru$variable))))
TB.Stru <- TB.Stru[TB.Stru$value != '',]

head(TB.Stru)
ggplot(TB.Stru, aes(x = human, y = variable, fill = value)) +
  geom_tile() +
  labs(x = 'Human Residue', y = 'Secondary Structure')
