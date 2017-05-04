################################################################################
# This script performs GSEA analysis on Expression matrix from Leukemia patients
# and outputs the enrichment result for ALL and AML cancer samples.
# Author: Suhas Vasaikar (suhasibb@gmail.com)

#To run script:
#1. Transfer input data into folder "dataset" (here leukemia.txt, pathways.txt, 
#   interestGenelist.txt)
#2. Create folder "output" for output result
#3. To run GSEA on local machine download GSEA tool from Broad (http://software.broadinstitute.org/gsea/downloads.jsp).
#4. Transfer gsea.jar into "jar" folder
################################################################################

################################################################################
#Step-1 : Input file validation
################################################################################
# Define path of directory
    PATH <- getwd(); #Change path if the GSEA and other files do not in the same folder
    timeStamp <- format(Sys.time(), "%m%d%Y_%H%M%S"); #Timestamp

# Change the input Expression matrix (leukemia.txt file) to GSEA comoatible format
# If matrix file is in the GSEA compatible forat ignore this step.
    temp <- read.table(paste(PATH, "/dataset/leukemia.txt",sep=""), check.names = F, header=T, sep="\t");
    NAME <- as.character(temp[,1]);
    DESCRIPTION <- 'na';
    data <- cbind(NAME, DESCRIPTION, temp[,2:ncol(temp)])
    write.table(data,"dataset/input_leukemiaExpression_forgsea.txt", row.names=F, col.names = T, quote=F, sep="\t");

# Create phenotype label file in GSEA comoatible format
# If phenotype file is in the GSEA compatible forat ignore this step.
    x1 <- dim(temp);                #temp variable, dimensions
    x1[2] <- x1[2]-1;               #number of actual labels
    x2 <- colnames(temp);           #temp variable, column names
    x3 <- unique(x2[2:length(x2)])  #temp variable, unique groups

    str1 <- paste(x1[2],"2","1\n",collapse = "\t");
    str2 <- paste(x2[2:length(x2)],collapse = "\t");
    str3 <- paste(x3,collapse = "\t");
    
    phenotypeFile <- paste(str1,'#',str3,'\n',str2, sep="")
    write.table(phenotypeFile,"dataset/input_leukemiaPhenotype_forgsea.cls", row.names=F, col.names = F, quote=F, sep="\t");

# Create Geneset file in GSEA comoatible format
# If geneset file is in the GSEA compatible forat ignore this step.
    #Open file "dataset/pathways.txt" and save it as "dataset/input_genesetPathwaysymbols.gmt"

#For batch query input can be defined here. Such as,
#    expFile <- "dataset/input_leukemiaExpression_forgsea.txt"    #expression matrix
#    phenFile <- "dataset/input_leukemiaPhenotype_forgsea.cls"    #phenotype label
#    genesetFile <- "dataset/input_genesetPathwaysymbols.gmt"     #pathway geneset

################################################################################

################################################################################
#Step-2 : Run GSEA on input data (leukemia.txt)
################################################################################
#Console
    #cmd <- "java -Xmx1024m -cp gsea2.jar xtools.gsea.Gsea -res Leukemia.gct -cls Leukemia.cls -gmx temp_pathwaysymbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 100 -scoring_scheme weighted -rpt_label gsea_result -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out GSEA_PHP-master/ -gui false";

#arguments for GSEA
    arg1 <- "-cp jar/gsea2.jar"                        #GSEA tool
    arg2 <- "-Xmx2048m"                                
    arg3 <- "xtools.gsea.Gsea"
    arg4 <- "-res dataset/input_leukemiaExpression_forgsea.txt"     #expression matrix
    arg5 <- "-cls dataset/input_leukemiaPhenotype_forgsea.cls"  #phenotype label
    arg6 <- "-gmx dataset/input_genesetPathwaysymbols.gmt"     #pathway geneset
    arg7 <- "-collapse false"
    arg8 <- "-mode Max_probe"
    arg9 <- "-norm meandiv"
    arg10 <- "-nperm 100"
    arg11 <- "-scoring_scheme weighted"
    arg12 <- "-rpt_label my_analysis"                  #Change report name
    arg13 <- "-include_only_symbols true"
    arg14 <- "-make_sets true"
    arg15 <- "-plot_top_x 20"                          #Top 20 result                
    arg16 <- "-rnd_seed 123"                           #Reproduce same result
    arg17 <- "-set_max 500"
    arg18 <- "-set_min 15"
    arg19 <- "-zip_report false"
    arg20 <- paste("-out gseaResult_",timeStamp,"/",sep="") #Output directry name
    arg21 <- "-gui false"
    cmd <- paste("java", arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20, arg21, sep=" ");
    
    #Run programe on console
    system(cmd);

    #Find Output folder
    folderName <- paste(PATH,"/gseaResult_",timeStamp,sep="");
    file.list <- as.character(list.files(folderName));
    fileNum <- gsub("my_analysis.Gsea.","",file.list[1]);

    #Read from Output folder and write into output folder
    ALL <- read.csv(paste(folderName,'/',file.list[1],"/gsea_report_for_ALL_",fileNum,".xls",sep=""), sep="\t");
    ALL <- ALL[order(ALL$NES, decreasing = T),];
    write.table(ALL,"output/gsea_pos_ALL.csv", row.names=F, col.names = T, sep=",");

    AML <- read.csv(paste(folderName,'/',file.list[1],"/gsea_report_for_AML_",fileNum,".xls",sep=""), sep="\t");
    AML <- AML[order(AML$NES, decreasing = F),];
    write.table(AML,"output/gsea_pos_AML.csv", row.names=F, col.names = T, sep=",");
################################################################################

################################################################################
#Step-3 : Find Gene function
################################################################################
    if (!require("qusage")) install.packages("qusage");
    if (!require("plyr")) install.packages("plyr");
    library("qusage");
    library("plyr");
    
    #Read geneset data
    geneSets = read.gmt("dataset/input_genesetPathwaysymbols.gmt");  #pathway geneset
    res_id <- ldply (geneSets, data.frame);
    colnames(res_id) <- c("Pathway", "Gene");

    #Read interest gene list and find the function mentioned in geneset
    interestGenelist <- read.csv("dataset/interestGenelist.txt", header=F)
    interestGene_pathwayres <- apply(interestGenelist,1,function(x) {
        gene <- as.character(x);
        temp <- res_id[res_id$Gene==gene,];
        if(nrow(temp)!=0) {
            res <- c(paste(as.character(temp[,1]),collapse =", "), nrow(temp));
        } else {
            res <- c('',0);
        }
        return(res);
    });
    interestGene_pathwayres <- data.frame(t(interestGene_pathwayres));
    interestGene_pathwayres  <- cbind(interestGenelist, interestGene_pathwayres);
    colnames(interestGene_pathwayres) <- c("Gene", "Pathway", "Count");
    write.table(interestGene_pathwayres,"output/interestGene_pathwayResult.csv", row.names=F, col.names = , sep=",");

    #All results are in output folder
    #output/gsea_pos_ALL.csv
    #output/gsea_pos_AML.csv
    #output/interestGene_pathwayResult.csv
#END
################################################################################
