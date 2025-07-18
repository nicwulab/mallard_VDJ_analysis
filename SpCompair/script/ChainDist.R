library(ggplot2)

Files <- list.files("result/ChainDist")

All = list()
for(F_id in Files){
  file.loc <- paste("result/ChainDist/", F_id, sep="")
  tb <- read.csv(file.loc)
  # extract the first letter form the column X
  SEQ = paste(substr(tb$X, 1, 1), collapse = "")
  # find the position of the first occurrence of "GGRPRQ"
  pos <- regexpr("GGRPRQ", SEQ) 
  # Match the position with the columns
  STR = c()
  for(i in c(pos[1]:(pos[1]+6))) {
    STR = c(STR, colnames(tb)[tb[i,] == min(tb[i, -1])])
  }
  All[[F_id]] <- STR
}


TB <- as.data.frame(t(as.data.frame(All)))

TBF <- TB[grep("distG", row.names(TB)),]
TBE <- TB[grep("distG", row.names(TB), invert = T),]


tb <- read.csv("result/ChainDist/fold_rag1_rag2_human_model_0_dist.csv")
# extract the first letter form the column X
SEQ = paste(substr(tb$X, 1, 1), collapse = "")
# find the position of the first occurrence of "GGRPRQ"
pos <- regexpr("GGRPRQ", SEQ) 

for(i in c(pos[1]:(pos[1]+6))) {
  print(colnames(tb)[tb[i,] == min(tb[i, -1])])
}


