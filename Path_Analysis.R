library(org.Hs.eg.db)
library(pathfindR)

Gene.symbol=row.names(resOrdered_padj)
logFC=resOrdered_padj[,2]
adj.P.Val=resOrdered_padj[,6]

input_df=data.frame(Gene.symbol, logFC,adj.P.Val)
head(input_df)
dim(input_df)

output_df <- run_pathfindR(input_df, gene_sets = "Reactome") #is possivble change gene_set by GO_all, etc
dim(output_df)


enrichment_chart(result_df = output_df, top_terms = 15) #PLOT

term_gene_graph(result_df = output_df, use_description = TRUE, num_terms = 14) #PLOT
