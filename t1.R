

print('hey')
colnames(tot)[11]<-'phen'
tot<-tot[which(tot$lod>8.38),]
