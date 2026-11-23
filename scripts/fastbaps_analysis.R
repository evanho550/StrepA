# =============================================================================
# FastBAPS Analysis — HPC version
# =============================================================================


library(fastbaps)
library(ggtree)
library(ggtreeExtra)
library(phytools)
library(ggplot2)
library(ape)
library(dplyr)

# =============================================================================
# File paths
# =============================================================================

fasta.file.name  <- file.path(Sys.getenv("PANAROO_RESULTS_DIR"), "ref_guided_core_gene_alignment_dedup.aln")
tree.file.name   <- file.path(Sys.getenv("RAXML_RESULTS_DIR"), "core_genome_search.raxml.bestTree")
output.dir       <- Sys.getenv("PANAROO_RESULTS_DIR")

# =============================================================================
# Step 1: Import alignment and optimise prior
# =============================================================================

message("Importing alignment...")
sparse.data <- import_fasta_sparse_nt(fasta.file.name)

message("Optimising prior (this may take >10 minutes)...")
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")

# =============================================================================
# Step 2: Cluster with FastBAPS
# =============================================================================

message("Running FastBAPS clustering...")
baps.hc        <- fast_baps(sparse.data)
best.partition <- best_baps_partition(sparse.data, baps.hc)

names(best.partition) <- colnames(sparse.data$snp.matrix)

# =============================================================================
# Step 3: Export cluster assignments to CSV
# =============================================================================

best.partition.df <- data.frame(
  isolate = names(best.partition),
  cluster = best.partition,
  row.names = NULL
)

print(head(best.partition.df))
print(table(best.partition.df$cluster))

write.csv(best.partition.df,
          file = file.path(output.dir, "fastbaps_clusters.csv"),
          row.names = FALSE)
message("Cluster assignments written to fastbaps_clusters.csv")

# =============================================================================
# Step 4: Load tree and verify sample ID matching
# =============================================================================

message("Loading tree and checking sample IDs...")
iqtree <- ape::read.tree(tree.file.name)
ids    <- colnames(sparse.data$snp.matrix)

if (!all(ids %in% iqtree$tip.label)) {
  stop("Sample IDs from alignment do not match tree tip labels. Check for naming inconsistencies.")
}

plot.df <- data.frame(
  id       = ids,
  fastbaps = as.factor(best.partition[ids]),
  stringsAsFactors = FALSE
)

# =============================================================================
# Step 5a: Visualisation — tile track alongside tree (facet_plot)
# =============================================================================

message("Generating tile track plot...")
gg <- ggtree(iqtree)

f2 <- facet_plot(
  gg,
  panel   = "FastBAPS Cluster",
  data    = plot.df,
  geom    = geom_tile,
  mapping = aes(x = 1, fill = fastbaps),
  color   = NA
) + scale_fill_viridis_d(name = "FastBAPS Cluster")

ggsave(file.path(output.dir, "fastbaps_tree_tile.pdf"),
       plot = f2, width = 12, height = 16)
message("Tile plot saved to fastbaps_tree_tile.pdf")

# =============================================================================
# Step 5b: Visualisation — coloured tip points on tree (geom_tippoint)
# =============================================================================

message("Generating tip point plot...")
tip_labels          <- iqtree$tip.label
cluster_assignments <- best.partition[tip_labels]

tippoint.df <- data.frame(
  id       = tip_labels,
  fastbaps = as.factor(cluster_assignments),
  stringsAsFactors = FALSE
)

p_tips <- ggtree(iqtree) %<+% tippoint.df +
  geom_tippoint(aes(color = fastbaps), size = 2) +
  scale_color_viridis_d(name = "FastBAPS Cluster") +
  theme_tree2()

ggsave(file.path(output.dir, "fastbaps_tree_tippoints.pdf"),
       plot = p_tips, width = 12, height = 16)
message("Tip point plot saved to fastbaps_tree_tippoints.pdf")

message("Analysis complete.")
