### Part 2 analyzes each guide in the screen for DEGs. 
# For these initial comparisons, the control group contains only non-targeting guides. 

# Define a list of non-targeting control guides. 
Control_guides <- c('non-targeting_00000','non-targeting_00001','non-targeting_00002','non-targeting_00003',
                    'non-targeting_00004','non-targeting_00005','non-targeting_00006','non-targeting_00007',
                    'non-targeting_00008','non-targeting_00009','non-targeting_00010','non-targeting_00011',
                    'non-targeting_00012','non-targeting_00013','non-targeting_00014','non-targeting_00015',
                    'non-targeting_00016')
All_guide_names <- rownames(CRISPR_guide_data)
# Non-control guides are the guides within All guides that are not considered Control as above
Non_control_guides <- setdiff(All_guide_names,Control_guides)
metadata_with_CRISPR$Control_guide_status <- 'Not_Control_Only'
# Label cells according to control guide status
for(cell in rownames(metadata_with_CRISPR)){
  if(any(metadata_with_CRISPR[cell,Control_guides] =='Contains')){
    if(!any(metadata_with_CRISPR[cell,Non_control_guides] == 'Contains')){
      metadata_with_CRISPR[cell,"Control_guide_status"] <- 'Control_Only'
    }
  }
}
# Turn this data into a factor, determine how many cells are control only 
metadata_with_CRISPR$Control_guide_status <- factor(metadata_with_CRISPR$Control_guide_status)
table(metadata_with_CRISPR$Control_guide_status) # 1021 with control only, 64855 else

### Next: Analyze DEGs on a per-guide basis. Take C3 enhancer locus 1 for example. 
# C3 sig e1 
# List of guides you want to iterate over
guides <- c('e1_g1','e1_g2','e1_g3','e1_g4','e1_g5') # Add all your guides here

# Iterate through each guide
for(guide in guides){
  
  # Dynamically create the column name
  guide_column <- paste0('C3_sig_', guide)  # E.g., 'C3_sig_e2_g1'
  guide_status_column <- paste0(guide_column, '_status')  # E.g., 'C3_sig_e2_g1_status'
  
  # Initialize the status column for each guide
  metadata_with_CRISPR[[guide_status_column]] <- paste0('No_', guide_column)
  
  # Update the status column based on the guide's condition
  for(cell in rownames(metadata_with_CRISPR)){
    if(any(metadata_with_CRISPR[cell, guide_column] == 'Contains')){
      metadata_with_CRISPR[cell, guide_status_column] <- guide_column
    }
  }
  
  # If control status is 'Control_Only', update the status to 'Control'
  metadata_with_CRISPR[metadata_with_CRISPR$Control_guide_status == 'Control_Only', guide_status_column] <- 'Control'
  
  # Convert the status column to a factor and print its levels
  metadata_with_CRISPR[[guide_status_column]] <- factor(metadata_with_CRISPR[[guide_status_column]])
  print(levels(metadata_with_CRISPR[[guide_status_column]]))
  
  # Update MPS_filtered with the new status column
  MPS_filtered@meta.data <- metadata_with_CRISPR
  MPS_filtered <- SetIdent(MPS_filtered, value = guide_status_column)
  
  # Find and store DEGs for the current guide vs control
  DEGs <- FindMarkers(MPS_filtered,
                      ident.1 = guide_column,
                      ident.2 = 'Control')
  # Write csvs of DEGs 
  write.csv(DEGs, file = paste0('mar5_guides_vs_NTC/',guide_column,'all_DEGs'))
}

# Further examination indicated the need for Benjamini-Hochberg p-value adjustments. 
#C3 Enhancer 1 Guide 1
C3_e1_g1_all_DEGs = read.csv('mar5_guides_vs_NTC/C3_e1_g1_all_DEGs.csv', row.names = 1)
p_values_C3_e1_g1 <- C3_e1_g1_all_DEGs$p_val # Original p-values
p_adjusted_C3_e1_g1 <- p.adjust(p_values_C3_e1_g1, method = "BH")
C3_e1_g1_all_DEGs$p_val_bh <- p_adjusted_C3_e1_g1
# Save the DEGs with the new BH p values
write.csv(C3_e1_g1_all_DEGs, file = paste0('mar25_bh_fc_filtered/C3_e1_g1_all_w_bh.csv'))
# Find DEGs with a p-value of less than 0.05
bh_sig_DEGs_C3_e1_g1 <- C3_e1_g1_all_DEGs[C3_e1_g1_all_DEGs$p_val_bh < 0.05, ]
# Find DEGs with at least a 15% fold change
bh_and_fc_DEGs_C3_e1_g1 <- bh_sig_DEGs_C3_e1_g1[abs(bh_sig_DEGs_C3_e1_g1$avg_log2FC) > 0.20163,]
write.csv(bh_and_fc_DEGs_C3_e1_g1, file = paste0('mar25_bh_fc_filtered/C3_e1_g1_bh_fc_filtered.csv'))