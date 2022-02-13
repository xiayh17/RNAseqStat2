devtools::load_all()
devtools::document()

# treatInfo(cutFC = NULL,label = c("Down", "Stable", "Up"),
#            label_ns = "Stable")
#
# rowSums(counts > 2) >= ncol(counts)

library(airway)

data("airway")

row_counts <- as.data.frame(assay(airway))

group_list <- as.character(colData(airway)$dex)

data_i <- Create_DEGContainer(expMatrix = row_counts,
                              groupInfo = group_list,
                              caseGroup = "trt",
                              idType = "ENSEMBL")
