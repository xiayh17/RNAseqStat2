test_that("test deg_DESeq2", {
  expect_message(deg_DESeq2(counts_input,group_list,
                            case_group = "T", control_group = "C", qc = TRUE,
                            x = "log2FoldChange", y = "pvalue",
                            dir = tempdir(), prefix = "2-DEG_DEseq2"),"were store in")
})
