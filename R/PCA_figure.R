
library(readxl)
install.packages("FactoMineR")
library("FactoMineR")


a_mat <- read_excel("CC/Supplementary table 4- Frequency Matrix First Set.xlsx")
missing <- read_excel("CC/Supplementary table 3- Missing Strains,first cohort.xlsx")

#find founder strains
index_founder <- which(sapply(a_mat$Strain,function(strain){
  length(grep("CC",strain))==0
}))

#fix founder B cell abundance
a_mat$Bcells[index_founder] <- a_mat$EarlyB[index_founder]+a_mat$LateB[index_founder]
write.csv(a_mat,file="CC/updated_abundance.csv")

#find non missing strains for PCA
rownames_missing <- missing$Strain
missing <- missing %>% select(-Strain)
rownames(missing)<-rownames_missing

non_missing <- rownames(missing)[which(rowSums(missing ) == ncol(missing))]

#subset the table to include non missing strains only 
a_mat <- a_mat[a_mat$Strain %in% non_missing,]
a_mat$idx <- paste0(a_mat$Strain,"_",1:nrow(a_mat))
rownames(a_mat)<- a_mat$idx
supp_table <- a_mat %>% select(Strain,idx)
a_mat <- a_mat %>% select(-c(Strain,idx))
rownames(a_mat)<-supp_table$idx

res.pca = PCA(a_mat, scale.unit = TRUE, ncp=2, graph=T)
pca_coord<- res.pca$ind$coord

res <- as.matrix(dist(pca_coord))

res <- reshape2::melt(res)
res$Var1<- as.character(res$Var1)
res$Var2<- as.character(res$Var2)

res$Var1<-sapply(res$Var1,function(x){strsplit(x,"_", fixed=TRUE)[[1]][1]})
res$Var2<-sapply(res$Var2,function(x){strsplit(x,"_", fixed=TRUE)[[1]][1]})

res <- res[-which(res$value==0),]

res$type <- "none"
for (r in 1:nrow(res)){
  
  if(res[r,"Var1"]==res[r,"Var2"]){
    res[r,"type"] <- "pairs"
  }
  if(grepl("CC",res[r,"Var1"]) && grepl("CC",res[r,"Var2"]) && res[r,"type"]=="none"){
    res[r,"type"] <- "CC"
  }
  if(!grepl("CC",res[r,"Var1"]) && !grepl("CC",res[r,"Var2"]) && res[r,"type"]=="none"){
    res[r,"type"] <- "founders"
  }
}

res <- res[-which(res$type %in% "none"),]
write.csv(res,"CC/distances.csv")
pca_coord <- as.data.frame(pca_coord)
pca_coord$"Strain"<-sapply(rownames(pca_coord),function(x){strsplit(x,"_", fixed=TRUE)[[1]][1]})

write.csv(pca_coord,"CC/pca_coord.csv")

