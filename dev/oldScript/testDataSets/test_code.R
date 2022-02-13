devtools::load_all()
plot <- enhance_volcano(deg_data_g,x = x, y = y, label = label, cut_FC = 3, cut_FDR = 0.001)
# plot
# #
# if (any(!is.null(genes_list),!is.null(highlight))) {
# plot <- show_genes(deg_data_g, x = x, y = y, genes_list, highlight,plot,size,expand)
# }
# any(!is.null(genes_list),!is.null(highlight))
plot +
  labs(colour = "Group",x = "log2FoldChange", y = "-log10(PValue)",
       subtitle = paste("hinihaoya \n niyehaoya"))
ggsave(plot,filename = "_volcano.pdf", width = 1600,height = 1600,units = "px",limitsize = FALSE)
