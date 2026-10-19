### Part 3 expands on the analyses in part 2 by defining a set of "pooled good guides," 
# or guides that have a significant impact on target gene expression. We also define 
# "no cis target" (NCT) as a group of guides, inclusive of non-targeting controls, for which
# we were unable to define a significant knockdown of the local cis target. 

# As an example, C3:
guides_no_cis_target <- c()

# Set your folder path
folder_path <- "apr11_bh_guides"

# List all files starting with the gene and ending with "bh.csv"
files <- list.files(folder_path, pattern = "^C3.*bh\\.csv$", full.names = TRUE)

# Loop through each file and check the first column
for (file in files) {
  data <- read.csv(file, stringsAsFactors = FALSE)
  
  # Extract first column
  first_col <- data[[1]]
  
  # Check for exact "gene", if not present, add guide to no cis target list
  if (!any(grepl("\\bC3\\b", first_col))) {
    base_name <- basename(file)
    name_without_suffix <- sub("_all_DEGs_w_bh\\.csv$", "", base_name)
    guides_no_cis_target <- c(guides_no_cis_target, name_without_suffix)
  }
}

# Repeat for all guides. 

# Non-control guides are the guides within All guides that are not considered Control as above
guides_w_cis_target <- setdiff(All_guide_names,guides_no_cis_target)
metadata_with_CRISPR$NCT_status <- 'Not_NCT_only'
# Label cells according to control guide status
for(cell in rownames(metadata_with_CRISPR)){
  if(any(metadata_with_CRISPR[cell,guides_no_cis_target] =='Contains')){
    if(!any(metadata_with_CRISPR[cell,guides_w_cis_target] == 'Contains')){
      metadata_with_CRISPR[cell,"NCT_status"] <- 'NCT_only'
    }
  }
}
# Turn this data into a factor
metadata_with_CRISPR$NCT_status <- factor(metadata_with_CRISPR$NCT_status)
table(metadata_with_CRISPR$NCT_status) # 9091 with NCT only, 56815 else

### C3 pooled good guides vs nocis

# Define a list of "good guides" that significantly impacted C3. 
C3_good_guides <- c('C3_sig_e2_g1','C3_sig_e2_g2','C3_sig_e2_g3','C3_sig_e2_g6','C3_sig_e2_g8',
                    'C3_tss_g1','C3_tss_g2')

# C3 Good Guides
metadata_with_CRISPR$C3_target_status <- 'C3_not_targeted'
for(cell in rownames(metadata_with_CRISPR)){
  if(any(metadata_with_CRISPR[cell,C3_good_guides] =='Contains')){
    metadata_with_CRISPR[cell,'C3_target_status'] <- 'C3_targeted'
  }
}

# Create a new metadata column that will track C3 good guides AND no cis target guides. 
metadata_with_CRISPR[['C3_targeted_vs_nocis']] <- 'C3_not_targeted'
metadata_with_CRISPR[metadata_with_CRISPR$NCT_status == 'NCT_only', 'C3_targeted_vs_nocis'] <- 'NCT_only'
metadata_with_CRISPR[metadata_with_CRISPR$C3_target_status == 'C3_targeted', 'C3_targeted_vs_nocis'] <- 'C3_targeted'
metadata_with_CRISPR$C3_targeted_vs_nocis <- factor(metadata_with_CRISPR$C3_targeted_vs_nocis)
MPS_filtered@meta.data <- metadata_with_CRISPR
MPS_filtered <- SetIdent(MPS_filtered, value = 'C3_targeted_vs_nocis')
C3_targeted_vs_nocis_DEGs <- FindMarkers(MPS_filtered,
                                         ident.1 = 'C3_targeted', ident.2 = 'NCT_only')
write.csv(C3_targeted_vs_nocis_DEGs, file = paste0('apr23_guides_vs_nocistarget/C3_targeted_vs_nocis_all_DEGs.csv'))
p_values_orig <- C3_targeted_vs_nocis_DEGs$p_val # Original p-values
p_adjusted <- p.adjust(p_values_orig, method = "BH")
C3_targeted_vs_nocis_DEGs$p_val_bh <- p_adjusted
# Find DEGs with a p-value of less than 0.05
bh_sig_cis_targets <- C3_targeted_vs_nocis_DEGs[C3_targeted_vs_nocis_DEGs$p_val_bh < 0.05, ]
# Find DEGs with at least a 15% fold change
bh_and_fc_cis_targets <- bh_sig_cis_targets[abs(bh_sig_cis_targets$avg_log2FC) > 0.20163,]
# Filter by pct expression
cis_targets_filtered_pct1 <- bh_and_fc_cis_targets[bh_and_fc_cis_targets$pct.1 > 0.1,]
cis_targets_filtered_pct2 <- cis_targets_filtered_pct1[cis_targets_filtered_pct1$pct.2 > 0.1,]
write.csv(cis_targets_filtered_pct2, file = paste0('apr23_guides_vs_nocistarget/C3_targeted_vs_nocis_all_DEGs_filtered.csv'))
