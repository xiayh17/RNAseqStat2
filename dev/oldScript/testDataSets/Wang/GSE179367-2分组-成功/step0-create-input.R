rm(list = ls())
counts_input = read.csv(file = "testDataSets/Wang/GSE179367-2分组-成功/GSE179367_gene_count.real.txt",
                        header = T,sep = "\t")
counts_input = fread(file = "testDataSets/Wang/GSE179367-2分组-成功/GSE179367_gene_count.real.txt",
                        header = T)
counts_input[1:4,1:4]



kp = !duplicated(counts_input$X)
table(kp)
counts_input=counts_input[kp,]
rownames(counts_input)=counts_input$X
counts_input=counts_input[,-1]
counts_input=floor(2^(counts_input+1))

group_list = rep(c("Ctrl","MMP9"),each =3);group_list

save(counts_input,group_list,file = 'input.Rdata')



