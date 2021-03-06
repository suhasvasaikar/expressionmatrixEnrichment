# ExpressionmatrixEnrichment
<b>Gene Set Enrichment based metabolic pathway analysis in Leukemia</b>

</br>This script performs GSEA analysis on Expression matrix from Leukemia patients
</br>and outputs the enrichment result for ALL and AML cancer samples.
</br> Author: Suhas Vasaikar (suhasibb@gmail.com)

</br><b>TASK</b>:
</br>The function calculate the gene set enrichment statistics for the list of gene sets from Leukemia expression matrix with the corresponding normalized enrichment scores and pvalues. Results should be sorted by the normalized enrichment score.
</br>Identidy the functions for the user based gene list.

</br><b>To run script:</b>
</br>1. Transfer input data into folder "dataset" (here leukemia.txt, pathways.txt, interestGenelist.txt)
</br>User can have their expression matrix in .gct, phenotype in .cls, and geneset in .gmt format (http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html). 
</br>2. Create folder "output" for output result
</br>3. To run GSEA on local machine download GSEA tool from Broad (http://software.broadinstitute.org/gsea/downloads.jsp).
</br>4. Transfer gsea.jar into "jar" folder
</br>5. To run the file, open Rstudio/R and setwd() to the folder. Run each step "click run" (prefferred for single dataset) or
</br>open terminal. set directory "expressionmatrixEnrichment-master/" and type
</br>$ Rscript gsea.R

</br><b>Output:</b>
</br>All results are in output folder
</br>output/gsea_pos_ALL.csv
</br>output/gsea_pos_AML.csv
</br>output/interestGene_pathwayResult.csv

Enjoy.
