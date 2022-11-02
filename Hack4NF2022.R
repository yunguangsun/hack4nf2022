install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))

library(devtools)
install_github("Bin-Chen-Lab/octad.db")
install_github("Bin-Chen-Lab/octad")

library(synapser)
library(tidyverse)
library(edgeR)
library(octad)

synLogin(email="Your email", password="Your password")
set.seed('99999')

#read the normalized PN primary tissue microarray data from synapse, original data from GSE14038
GSE14038 <- synGet("syn5950004")
GEO<-read.table(GSE14038$path, header=TRUE,row.names=1,stringsAsFactors = FALSE,sep="\t")

GSE14038.anno<- synGet("syn5950620")
GEO.anno<-read.table(GSE14038.anno$path, header=TRUE,row.names=1,stringsAsFactors = FALSE,sep="\t")

#Generate a expression set
library(Biobase)
pData<-GEO.anno[order(rownames(GEO.anno)),]
pDataPNF<-new("AnnotatedDataFrame",data=pData)
GEO.eSet<-ExpressionSet(assayData=as.matrix(GEO),phenoData=pDataPNF,annotation="GSE14038")


# collect plexiform neurofibromas and Schwann cell transcriptome from the data set.
CollectedSample<-GEO.eSet[,GEO.eSet$cellType %in% c("pNF","NHSC")]


#Generate the differentially expressed gene between Plexifiorm neurofibromas and normal human Schwann cells.
#We hypothesize that decrease the highly expressed genes in tumors will inhibit the pN growth
library(edgeR)
y <- DGEList(counts=exprs(CollectedSample))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, ,keep.lib.sizes=FALSE]
ExpMatrix.trim <- y


group<-factor(CollectedSample$cellType)
y <- DGEList(counts=exprs(CollectedSample),group=group)

y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

et<-exactTest(y,pair=c("NHSC","pNF"))

Result.PNF<-topTags(et, n=20000, adjust.method="BH", sort.by="p.value")
pNF.DEGs<-data.frame(Result.PNF[(Result.PNF$table$FDR<0.05)&((Result.PNF$table$logFC>1)|(Result.PNF$table$logFC<(-1))),])

pNF.DEGs$Symbol<-rownames(pNF.DEGs)

#Format the result file for downstream analysis using the octad package. sRGES.pNF has the predicted drug targets from the LINCS database
colnames(pNF.DEGs)<-c("log2FoldChange","logCPM", "PValue","FDR","Symbol")
sRGES.pNF=runsRGES(pNF.DEGs,max_gene_size=100,permutations=1000,output=TRUE, outputFolder = "Your file path" )

#Add a new column "std_name" to the sRGES.pNF for the drug names
sRGES.pNF$std_name<-toupper(sRGES.pNF$pert_iname)



#Please assign a folder in your system to save the result.
#This process is to evaluate the chembl_targets enrichment of the predicted drugs
resultsDrug.pNF<-octadDrugEnrichment(sRGES = sRGES.pNF, target_type='chembl_targets',outputFolder = "Your file path")



#Download the drug screening results
drug_data <- synGet("syn20700260")$path %>% read.csv()

#Seven plexiform neurofibroma related cells were used in the drug screen.
CellLines <- c("ipNF05.5 (single clone)", "ipNF06.2A", "ipNF95.11b C/T", "ipnNF95.11C", "ipNF95.6", "ipNF05.5 (mixed clone)", "ipNF95.11b C")

#select the drug that has response data from PN cell lines.
drug_data_filt_1 <- drug_data %>%
  filter(response_type == "AUC_Simpson") %>%
  filter(model_name %in% CellLines) %>%
  group_by(drug_screen_id) %>%
  filter(n() == 1) %>%
  ungroup()


#narrow down the drugs have response data from at least three PN cell lines and generate a new column named "median_response".
drug_data_filt <- drug_data_filt_1 %>%
  group_by(DT_explorer_internal_id) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  dplyr::select(DT_explorer_internal_id, response) %>%
  group_by(DT_explorer_internal_id) %>%
  summarize('median_response' = median(response))%>%
  ungroup()


drug_data_filt_2 <- drug_data_filt_1 %>%
  dplyr::select(drug_name,DT_explorer_internal_id, response) %>%
  group_by(drug_name) %>%
  summarize('median_response' = median(response))%>%
  ungroup()  

#Download the drug annotation file from synapse and narrow down the drug candidates based on the pchembl database annotation.
targets <- synGet("syn17091507")$path %>% readRDS() %>%
  filter(mean_pchembl > 6) %>%
  dplyr::select(internal_id, hugo_gene, std_name) %>%
  distinct()

#Create a index file to connect DT_explorer_internal_id and drug_name
AnnoDrug <- drug_data_filt_1%>%dplyr::select(DT_explorer_internal_id,drug_name)%>%
  group_by(DT_explorer_internal_id)%>%unique()


#Create a table has median_response, hugo_gene, and std_name
DrugNameRes.final<- drug_data_filt_2 %>% left_join(AnnoDrug, by= c("drug_name"="drug_name")) %>% left_join(targets, by= c("DT_explorer_internal_id"="internal_id"))%>%
  dplyr::select(drug_name,DT_explorer_internal_id,median_response,hugo_gene,std_name)                  



#We overlap the predicted drug candidates with drug used in the PN drug screen
Overlap.pNF.prediction<-DrugNameRes.final%>% left_join(sRGES.pNF, by= c("std_name"="std_name"))%>%
  dplyr::select(drug_name,,median_response,hugo_gene,std_name,pert_iname,sRGES)  

#sRGES score less than -0.2 are recommended by the octad authors for drugs can reverse the increased gene (Nature Protocols volume 16, pages728–753 (2021))
# We set an arbitrary cut off in the PN drug screen to select drugs that have median_response less than 60% among PN cell lines.

Overlap.pNF.prediction.filter<-Overlap.pNF.prediction%>%
  filter(median_response<60)%>%filter(sRGES<(-.2)) %>% arrange(median_response)

# The overlapped drugs are the top candidates that can reverse the increased genes in PN tumors and are also verified by the drug screen using the PN cell lines.
# We recommend prioritizing these candidates in further experiment verification.

setwd("Your file path")
write.table( Overlap.pNF.prediction.filter,"Overlap.pNF.prediction.filter.txt",sep="\t")


