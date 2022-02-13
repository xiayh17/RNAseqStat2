## code to prepare `counts_input` counts_input goes here

load("data-raw/counts_input.Rdata")
write.table(group_list,file = "inst/extdata/group_list.txt",row.names = F,
            col.names = F,quote = F)
write.table(counts_input,file = "inst/extdata/counts_input.tsv",quote = F)
filename = "inst/extdata/counts_input.tsv"
R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=FALSE, remove=TRUE, BFR.SIZE=1e+07)
group_list <- readLines("inst/extdata/group_list.txt")
filename = "C:/Users/xiayh17/Desktop/RNAseqStat_STANTARDWORKFLOW/counts_input.tsv"
R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=FALSE, remove=TRUE, BFR.SIZE=1e+07)
counts_input <- read.csv("inst/extdata/counts_input.tsv.gz",sep = "")

# system.file("extdata", "mydata.csv", package = "RNAseqStat")
