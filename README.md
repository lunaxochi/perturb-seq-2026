# The overall purpose of these scripts is to analyze Perturb-seq data from iMGLs.
# Part 1 takes us from feature-barcode matrices to a unified Seurat object. 
# It also filters for standard QC metrics and creates a metadata table with CRISPR guide information. 
# Part 2 analyzes each guide in the screen for DEGs. 
# For these initial comparisons, the control group contains only non-targeting guides. 
# Part 3 expands on the analyses in part 2 by defining a set of "pooled good guides," 
# or guides that have a significant impact on target gene expression. We also define 
# "no cis target" (NCT) as a group of guides, inclusive of non-targeting controls, for which
# we were unable to define a significant knockdown of a local cis target. 
