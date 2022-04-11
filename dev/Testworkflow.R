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
                              dataType = "Array",
                              caseGroup = "trt",
                              idType = "ENSEMBL")

data_o <- runALL(object = data_i,dir = "output_test")
# data_t <- runMSigDB(obj = data_o,dir = "output_test")

# runCheck(data_i)
#
# data_g <- runDEG(obj = data_i, parallel = T)
#
# data_h <- runHyper(obj = data_g)
#
# data_gse <- runGSEA(obj = data_h)
#
# data_gsva <- runMSigDB(data_gse)
