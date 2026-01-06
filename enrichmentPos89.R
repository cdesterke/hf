list.files()

df<-read.table("pamr190.csv",h=T)


library(dplyr)

df%>%filter(HF.score>0)->pos

vector<-pos$id


# Installer Bioconductor et gprofiler2 si nécessaire
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("gprofiler2")

# Charger le package
library(gprofiler2)

# Exemple : liste de gènes humains
genes <- c("TP53", "BRCA1", "EGFR", "MYC")

# Analyse d’enrichissement avec g:Profiler
res <- gost(query = vector, organism = "hsapiens")

# Afficher les premiers résultats
head(res$result)

# Visualisation des pathways enrichis
gostplot(res, capped = TRUE, interactive = TRUE)

# Exporter les résultats en tableau
library(tidyr)



df <- as.data.frame(res$result)
df_flat <- unnest(df, cols = everything())




write.csv(df_flat, "gprofiler_results.csv", row.names = FALSE)

library(pheatmap)

df <- as.data.frame(res$result)

# Exemple : matrice avec -log10(p_value)
mat <- matrix(-log10(df_flat$p_value), nrow = 1)
colnames(mat) <- df_flat$term_name



top_terms <- head(df_flat[order(df_flat$p_value), ], 10)

mat <- as.matrix(top_terms[, c("p_value")])
rownames(mat) <- top_terms$term_name

pheatmap(-log10(mat), cluster_cols = FALSE,
         main = "Top 10 termes enrichis")



## top 5 par database

library(dplyr)
library(pheatmap)

# On suppose que df_flat est déjà un data.frame avec colonnes :
# term_name, p_value, source (la base de données : GO:BP, KEGG, Reactome, etc.)

# Sélection du top5 par base de données
top_terms_by_db <- df_flat %>%
  group_by(source) %>%
  arrange(p_value, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

# Création de la matrice : -log10(p_value)
mat <- as.matrix(top_terms_by_db$p_value)
rownames(mat) <- paste(top_terms_by_db$source, top_terms_by_db$term_name, sep = "_")

# Transformation en matrice avec -log10(p_value)
mat <- matrix(-log10(top_terms_by_db$p_value),
              ncol = 1,
              dimnames = list(rownames(mat), "significance"))

# Heatmap
pheatmap(mat,
         cluster_cols = FALSE,
         cluster_rows = F)




df <- as.data.frame(res$result)
df_subset <- df_flat[, c("term_name", "source", "intersection")]

head(df_subset, 5)

library(pheatmap)
pheatmap(mat, cluster_cols = TRUE, main = "Enrichissement g:Profiler")


## dotplot cleveland

library(ggplot2)
library(dplyr)
library(stringr)

# Sélection du top10 par base de données
top_terms_by_db <- df_flat %>%
  group_by(source) %>%
  arrange(p_value, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(term_short = str_trunc(term_name, 40))   # raccourcir les titres

# Cleveland dotplot en 2 colonnes de 5
ggplot(top_terms_by_db,
       aes(x = -log10(p_value),
           y = reorder(term_short, -log10(p_value)),
           size = intersection_size,
           color = source)) +
  geom_point(alpha = 1) +
  facet_wrap(~source, scales = "free_y", ncol = 2) +   # 2 colonnes
  labs(x = "-log10(p-value)",
       y = "enriched functions",
       size = "Number of genes",
       color = "databse") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

## clusterprofiler

library(clusterProfiler)
library(org.Hs.eg.db)


ego <- enrichGO(gene= vector, OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH")

pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
                
                
  ego <- enrichGO(gene= row.names(pos), OrgDb= org.Mm.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",qvalueCutoff  = 0.05)              
                




dotplot(ego, showCategory=10)
library(DOSE)
#BiocManager::install("enrichplot")
library(enrichplot)

genelist<-split(pos$HF.score,pos$id)
genelist <- as.numeric(genelist)

cnetplot(ego,foldChange=genelist,colorEdge=TRUE,showCategory=10)

ego2 <- pairwise_termsim(ego)
treeplot(ego2)

upsetplot(ego)



# Extraire les gènes séparés par "/"
term_gene_long <- do.call(rbind, lapply(1:nrow(top10), function(i) {
  genes <- strsplit(top10$geneID[i], "/")[[1]]
  data.frame(
    term_id = top10$ID[i],
    term = top10$Description[i],
    p.adjust = top10$p.adjust[i],
    Count = top10$Count[i],
    gene = genes,
    stringsAsFactors = FALSE
  )
}))

head(term_gene_long)


# Aperçu
head(term_gene_long)





library(igraph)

# 1) Construction du graphe bipartite (nœuds: termes + gènes, arêtes: appartenance)
edges_bip <- term_gene_long %>%
  dplyr::select(term, gene)



g_bip <- graph_from_data_frame(edges_bip, directed = FALSE)

# 2) Marquage bipartite (type = TRUE pour termes, FALSE pour gènes)
V(g_bip)$type <- V(g_bip)$name %in% unique(term_gene_long$term)

# 3) Mesures de centralité (sur les gènes uniquement)
is_gene <- !V(g_bip)$type
deg <- degree(g_bip)
bet <- betweenness(g_bip, directed = FALSE, normalized = TRUE)
cls <- closeness(g_bip, normalized = TRUE)

centrality_bip_genes <- tibble(
  node = V(g_bip)$name,
  is_gene = is_gene,
  degree = deg,
  betweenness = bet,
  closeness = cls
) %>%
  filter(is_gene) %>%
  arrange(desc(degree), desc(betweenness))

# Top hubs (bipartite)
head(centrality_bip_genes, 15)





# Exports
write.csv(term_gene_long, "top10_term_gene_edges.csv", row.names = FALSE)
write.csv(centrality_bip_genes, "centrality_bipartite_genes.csv", row.names = FALSE)


# Petite visualisation des 30 gènes les plus connectés (bipartite)
library(ggplot2)

top_bip <- centrality_bip_genes %>% slice_head(n = 30)


library(ggplot2)




ggplot(top_bip, aes(x = reorder(node, degree), y = degree)) +
  # Lignes horizontales (de 0 jusqu'au point)
  geom_segment(aes(x = reorder(node, degree), xend = reorder(node, degree),
                   y = 0, yend = degree),
               color = "gray70") +
  # Points dimensionnés par betweenness et colorés par closeness
  geom_point(aes(size = betweenness, color = closeness), alpha = 1) +
  coord_flip() +
  labs(x = "Genes",
       y = "Degrees (GO::BP connections)",
       title = "Top hub genes (GO::BP) in heart failure",
       size = "Betweenness",
       color = "Closeness") +
  theme_classic(base_size=14) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold")
  ) +
  scale_color_gradient(low = "royalblue", high = "darkred")

