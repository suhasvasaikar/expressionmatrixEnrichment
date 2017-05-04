# expressionmatrixEnrichment
Gene Set Enrichment based metabolic pathway analysis in Leukemia

</br>This script performs GSEA analysis on Expression matrix from Leukemia patients
</br>and outputs the enrichment result for ALL and AML cancer samples.
</br> Author: Suhas Vasaikar (suhasibb@gmail.com)

</br>To run script:
</br>1. Transfer input data into folder "dataset" (here leukemia.txt, pathways.txt, interestGenelist.txt)
</br>User can have their expression matrix in .gct, phenotype in .cls, and geneset in .gmt format (http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html). 
</br>2. Create folder "output" for output result
</br>3. To run GSEA on local machine download GSEA tool from Broad (http://software.broadinstitute.org/gsea/downloads.jsp).
</br>4. Transfer gsea.jar into "jar" folder

</br><b>Output:</b>
</br>All results are in output folder
</br>output/gsea_pos_ALL.csv
</br>output/gsea_pos_AML.csv
</br>output/interestGene_pathwayResult.csv

Enjoy.
