rm(list = ls())
fs = list.files(path = "testDataSets/Wang/GSE184648-2分组3重复-ensembl的ID/GSE184648_RAW/",full.names = T);
fs

library(data.table)
mat =do.call(cbind,
             lapply(fs, function(x){
               fread(x,data.table = F)[,7]
             }))
gid = fread(fs[1],data.table = F)[,1]

mat[1:4,1:4]
dim(mat)


library(AnnoProbe)
ids = annoGene(gid,'ENSEMBL')
ids[1:4,1:4]
dim(ids)
length(unique(ids$SYMBOL))
gs = ids[match(gid,ids$ENSEMBL),1]
length(unique(gs))

kp= !duplicated(gs)
table(kp)
counts_input = mat[kp,]
rownames(counts_input) = gs[kp]
counts_input[1:4,1:4]

keep_feature <- rowSums (counts_input > 2) > 4
table(keep_feature)
counts_input = counts_input[keep_feature,]
counts_input[1:4,1:4]

group_list = rep(c("p66","p71"),each =3);group_list
save(counts_input,group_list,file = 'testDataSets/Wang/GSE184648-2分组3重复-ensembl的ID/input.Rdata')
