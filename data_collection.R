## DATA COLLECTION FOR EPathDash ##

# This scripts collects the relevant information from the compendium of RNA-seq datasets used 
# in the CF-Seq application. The resulting Rdata file is used in the pathway_activation_analysis.R
# script to create the database for the E.PathDash application. Use of the external gene ID 
# mapping tool from the Uniprot Consortium is required (https://www.uniprot.org/id-mapping). See section
# "UniProt Translation" in this script for required file names for results returned by the tool.

# libraries
library(stringr)
library(readxl)
library(tidyverse)

# datadir defines the path to the Rdata file generated for the CF-Seq application
# instructions can be found here: https://github.com/DartCF/cf-seq
datadir="change_me"
# outdir defines location to write temporary files used for gene ID translation
outdir="change_me"
# transdir defines location for results from gene ID translation using the Uniprot Consortium translation tool
transdir="change_me"
# cfseq_dir defines location of UniProt annotations collected for CF-Seq (from https://github.com/DartCF/cf-seq)
cfseq_dir="change_me"
# new_datadir defines location of new data files for application 
new_datadir="change_me"

# data file from CF-Seq application (data objects: CFSeq_Data and GenePathway_Data)
load(paste(datadir,"Data.Rdata",sep=""))

# New Data Load -----------------------------------------------------------
# Before running this, run the pre_data_collection script to check the availability of gene IDS in the current
# version of the gene id dictionary.
#
# This code only handles new datasets of streptococcus, staphylococcus, and pseudomonas. Update the code to handle
# new species
{
  new_datafile<-"New_Data_File"
  load(paste0(new_datadir,new_datafile))
  
  species<-"Streptococcus_Species"
  new_studies<-names(New_EPathDash_Data[[species]])
  for(i in 1:length(new_studies)){
    study<-new_studies[i]
    CFSeq_Data[[species]][study]<-New_EPathDash_Data[[species]][i]
    print(CFSeq_Data[[species]][[study]][["Additional_Metadata"]]$GEO.Accession)
  }
  species<-"Staphylococcus_aureus"
  new_studies<-names(New_EPathDash_Data[[species]])
  for(i in 1:length(new_studies)){
    study<-new_studies[i]
    CFSeq_Data[[species]][study]<-New_EPathDash_Data[[species]][i]
    print(CFSeq_Data[[species]][[study]][["Additional_Metadata"]]$GEO.Accession)
  }
  species<-"Pseudomonas_aeruginosa"
  new_studies<-names(New_EPathDash_Data[[species]])
  for(i in 1:length(new_studies)){
    study<-new_studies[i]
    CFSeq_Data[[species]][study]<-New_EPathDash_Data[[species]][i]
    print(CFSeq_Data[[species]][[study]][["Additional_Metadata"]]$GEO.Accession)
  }
}
# Accessory Variables -----------------------------------------------------
# Study Information
species_list=c("Pseudomonas_aeruginosa","Bacteroides_Species","Staphylococcus_aureus","Streptococcus_Species")
pseudomonas_studies=list("GSE124385","GSE125646","GSE142448","GSE142464","GSE148597","GSE156995","GSE166602",
                         "GSE87213","GSE110445","GSE123356","GSE125646","GSE130190","GSE136111","GSE138731",
                         "GSE152480","GSE153067","GSE163234","GSE163248","GSE166986","GSE179150","GSE71880")
bacteroides_studies=list("GSE129572")
staph_studies=list("GSE122048","GSE125741","GSE130777","GSE132179","GSE148024","GSE89791")
strep_studies=list("GSE89964","GSE90021","GSE97218","GSE97357")

# New studies added from additional data load
new_pseudomonas_studies=names(New_EPathDash_Data[["Pseudomonas_aeruginosa"]])
new_staph_studies=names(New_EPathDash_Data[["Staphylococcus_aureus"]])
new_strep_studies=names(New_EPathDash_Data[["Streptococcus_Species"]])
pseudomonas_studies<-c(pseudomonas_studies,new_pseudomonas_studies)
staph_studies <- c(staph_studies,new_staph_studies)
strep_studies <- c(strep_studies,new_strep_studies)

queries=list(pseudomonas_studies,bacteroides_studies,staph_studies,strep_studies)
names(queries)=species_list

# Query data for relevant studies
App_Data<-list()
Raw_Count_Data<-list()
for (species in names(queries)){
  for (study in queries[[species]]){
    study_data<-list()
    study_data[["Additional_Metadata"]]<-CFSeq_Data[[species]][[study]][["Additional_Metadata"]]
    full_results_idx<-str_detect(names(CFSeq_Data[[species]][[study]]),"Full Results")
    study_data<-append(study_data,CFSeq_Data[[species]][[study]][full_results_idx])
    App_Data[[species]][[study]]<-study_data
    Raw_Count_Data[[study]][["Count_Table"]]<-CFSeq_Data[[species]][[study]][["Count_Table"]]
    Raw_Count_Data[[study]][["Design_Matrix"]] <- CFSeq_Data[[species]][[study]][["Design_Matrix"]]                           
  }
}

strains=list() # list of strains by study
for (species in names(App_Data)){
  strains=append(strains,lapply(App_Data[[species]],function(x){
    x[["Additional_Metadata"]][["Strain"]]
  }))
}

gene_ids=list() # list of gene identifiers used in expression dataframes, grouped by study 
for (species in names(App_Data)){
  for (study in names(App_Data[[species]])){
    comparison<-names(App_Data[[species]][[study]][2])[1]
    gene_ids[[study]]<-rownames(App_Data[[species]][[study]][[comparison]])
  }
}


# Data for UniProt Translation --------------------------------------------

#-- RS gene translation dictionaries
# Pseudomonas: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL23199&id=71102&db=GeoDb_blob151
#   Pseudomonas_RS_translate.txt
# Staph: https://aureowiki.med.uni-greifswald.de/Main_Page
#   GeneSpecificInformation_Newman.xlsx
#   GeneSpecificInformation_USA300_FPR3757.xlsx
# Strep: https://frontiersin.figshare.com/articles/dataset/Table_1_Manganese_Depletion_Leads_to_Multisystem_Changes_in_the_Transcriptome_of_the_Opportunistic_Pathogen_Streptococcus_sanguinis_XLSX/13192721/1
#   GSE150593_SK36_Translation.XLSX
{
  # file from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL23199&id=71102&db=GeoDb_blob151
  RS_gene_dict<-read.table(paste(getwd(),"Pseudomonas_RS_translate.txt",sep="/"),sep="\t",header = T)
  RS_gene_dict<-RS_gene_dict %>% distinct(old_locus_tag,new_locus_tag,.keep_all = TRUE)
  
  # dictionary for NWMN/USA300 gene translations (https://aureowiki.med.uni-greifswald.de/download_gene_specific_information)
  NWMN_RS_dict<-read_xlsx(paste(getwd(),"GeneSpecificInformation_Newman.xlsx",sep="/"), col_names = T)
  NWMN_UniProt_dict<-read.csv(paste0(cfseq_dir,"SA_Newman.csv"),header = T) %>% 
    select(UNIPROTKB,KEGG) %>% 
    mutate(GENENAME=str_split_i(KEGG,":",2))
  USA300_RS_dict<-read_xlsx(paste(getwd(),"GeneSpecificInformation_USA300_FPR3757.xlsx",sep="/"),col_names = T)
  USA300_UniProt_dict <- read.csv(paste0(cfseq_dir,"SA_USA300.csv"),header = T) %>% 
    select(UNIPROTKB, GENENAME)
  
  # dictionary for SK36 gene translations
  SK36_RS_dict<-read_xlsx(getwd(),"GSE150593_SK36_Translation.XLSX",sheet=2)
  SK36_UniProt_dict<-read.csv(paste0(cfseq_dir,"SS_SK36.csv"),header = T) %>% 
    select(UNIPROTKB,KEGG) %>% 
    mutate(GENENAME=str_split_i(KEGG,":",2)) %>% 
    select(-KEGG)
}

# GSE125646 
gse125646<-gene_ids[["GSE125646"]]

# GSE142448
write.csv(gene_ids[["GSE142448"]],file=paste(outdir,"gse142448.txt",sep=""),row.names = F,quote = F)

# GSE142464
gse142464_RS<-gene_ids[["GSE142464"]][str_detect(gene_ids[["GSE142464"]],"RS")==TRUE]
write.csv(gse142464_RS,file=paste(outdir,"gse142464_rs.txt",sep=""),row.names = F, quote = F)
gse142464_PA14<-RS_gene_dict[RS_gene_dict$old_locus_tag %in% gse142464_RS,]$new_locus_tag
gse142464_PA14<-gse142464_PA14[nzchar(gse142464_PA14)]
gse142464_gene_name<-gene_ids[["GSE142464"]][!(str_detect(gene_ids[["GSE142464"]],"PA14"))]

# GSE124385
gse124385_RS<-gene_ids[["GSE124385"]][str_detect(gene_ids[["GSE124385"]],"RS")==TRUE]
gse124385_PA14<-RS_gene_dict[RS_gene_dict$old_locus_tag %in% gse124385_RS,]$new_locus_tag
gse124385_PA14<-gse124385_PA14[nzchar(gse124385_PA14)]

# GSE148597
gse148597_PA14<-gene_ids[["GSE148597"]][str_detect(gene_ids[["GSE148597"]],"PA14")]
gse148597_gene_name<-gene_ids[["GSE148597"]][!(str_detect(gene_ids[["GSE148597"]],"PA14"))]

# GSE156995
gse156995<-gene_ids[["GSE156995"]][str_detect(gene_ids[["GSE156995"]],"PA14")]

# GSE16602
write.csv(gene_ids[["GSE166602"]],file=paste(outdir,"gse166602.txt",sep=""),row.names = F,quote = F)

# GSE87213
gse87213_PA14<-gene_ids[["GSE87213"]][str_detect(gene_ids[["GSE87213"]],"PA14")]
gse87213_EBG<-gene_ids[["GSE87213"]][str_detect(gene_ids[["GSE87213"]],"EBG")]

# GSE110445
write.csv(gene_ids[["GSE110445"]],file=paste(outdir,"gse110445.txt",sep=""),row.names = F,quote = F)

# GSE123356
write.csv(gene_ids[["GSE123356"]],file=paste(outdir,"gse123356.txt",sep=""),row.names = F,quote = F)

# GSE130190
gse130190_PA<-gene_ids[["GSE130190"]][str_detect(gene_ids[["GSE130190"]],"PA")]
gse130190_gene_name<-gene_ids[["GSE130190"]][!(str_detect(gene_ids[["GSE130190"]],"PA"))]

# GSE136111
gse136111_PA<-gene_ids[["GSE136111"]][str_detect(gene_ids[["GSE136111"]],"PA")]
gse136111_gene_name<-gene_ids[["GSE136111"]][!(str_detect(gene_ids[["GSE136111"]],"PA"))] #UTR_XXXX don't know how to map

# GSE138731
write.csv(gene_ids[["GSE138731"]],file=paste(outdir,"gse138731.txt",sep=""),row.names = F,quote = F)

# GSE152480 --> dropped
gene_name_idx<-which(str_detect(gene_ids[["GSE152480"]],"PA")==FALSE)
gse152480_PA<-lapply(gene_ids[["GSE152480"]],function(x){str_extract(x,"^CDS\\.PA[0-9]+")})
gse152480_PA<-unlist(lapply(gse152480_PA,function(x){str_extract(x,"PA[0-9]+")})) # don't omit NA for Uniprot data struct
gse152480_PA<-na.omit(unlist(lapply(gse152480_PA,function(x){str_extract(x,"PA[0-9]+")}))) # omit NA for Uniprot translation
gse152480_gene_name<-str_split(gene_ids[["GSE152480"]][gene_name_idx][2],"\\.")[[1]][4]

# GSE153067 --> dropped
gse153067_PA<-lapply(gene_ids[["GSE153067"]],function(x){str_extract(x,"PA[0-9]+\\.[0-9]*")})
gse153067_PA<-unlist(lapply(gse153067_PA,function(x){
  if (nchar(x) == 7){
    return(substr(x,1,6))
  }else{
    return(x)
  }
}))

# GSE163234 --> dropped, strain PAO1161 not available in Uniprot mapper

# GSE163248
gse163248_PA<-gene_ids[["GSE163248"]][str_detect(gene_ids[["GSE163248"]],"PA")]
gse163248_gene_name<-gene_ids[["GSE163248"]][!(str_detect(gene_ids[["GSE163248"]],"PA"))]

# GSE166986 --> dropped
gse166986<-read_excel(paste(getwd(),"/gse166986_supplemental_files/Table S5. RNA-seq comparisons.xlsx",sep=""))
gse166986<-gse166986[2:nrow(gse166986),c("CDS number","PA number")]

# GSE179150
gse179150_PA<-gene_ids[["GSE179150"]][str_detect(gene_ids[["GSE179150"]],"PA")]
gse179150_gene_name<-gene_ids[["GSE179150"]][!(str_detect(gene_ids[["GSE179150"]],"PA"))]

# GSE71880
gse71880_PA<-gene_ids[["GSE71880"]][str_detect(gene_ids[["GSE71880"]],"PA")]

# GSE185398
gse185398_RS<-gene_ids[["GSE185398"]][str_detect(gene_ids[["GSE185398"]],"RS")==TRUE]
gse185398_PA14<-RS_gene_dict[RS_gene_dict$old_locus_tag %in% gse185398_RS,]$new_locus_tag
gse185398_PA14<-gse185398_PA14[nzchar(gse185398_PA14)]

# Pseudomonas_aeruginosa
PAO1_Numbers<-union(union(union(union(union(union(union(union(union(gse71880_PA,gse179150_PA),
                                            gse166986$`PA number`),
                                      gse163248_PA),
                                gse153067_PA),
                          gse152480_PA),
                    gse136111_PA),
                  gse130190_PA),
                  gene_ids[["GSE123356"]]),
                  gene_ids[["GSE110445"]])
PAO1_gene_names<-union(union(union(union(gse130190_gene_name,gse136111_gene_name),
                                   gse152480_gene_name),
                             gse163248_gene_name),
                       gse179150_gene_name)
PA14_Numbers<-union(union(union(union(union(union(union(union(gse142464_PA14,gse148597_PA14),
                                 gse156995),
                           gene_ids[["GSE166602"]]),
                     gse87213_PA14),
                     gene_ids[["GSE142448"]]),
                     gse125646),
                     gse185398_PA14)
PA14_gene_names<-union(gse142464_gene_name,gse148597_gene_name)

write.csv(PAO1_Numbers,file=paste(outdir,"PA01_Numbers.txt",sep=""),row.names = F,quote = F)
write.csv(PAO1_gene_names,file=paste(outdir,"PAO1_gene_names.txt",sep=""),row.names = F,quote = F)
write.csv(PA14_Numbers,file=paste(outdir,"PA14_Numbers.txt",sep=""),row.names = F,quote = F)
write.csv(PA14_gene_names,file=paste(outdir,"PA14_gene_names.txt",sep=""),row.names = F,quote = F)

# Bacteroides_Species
gse129572_BT<-gene_ids[["GSE129572"]][str_detect(gene_ids[["GSE129572"]],"BT")]
write.csv(gse129572_BT,file = paste(outdir,"BT_Numbers.txt",sep=""),row.names = F,quote = F)

# Staphylococcus_aureus
gse122048_NWM<-gene_ids[["GSE122048"]][str_detect(gene_ids[["GSE122048"]],"NWMN")]
gse122048_gene_name<-gene_ids[["GSE122048"]][!(str_detect(gene_ids[["GSE122048"]],"NWMN"))]
gse125741_USA<-gene_ids[["GSE125741"]][str_detect(gene_ids[["GSE125741"]],"SAUSA")]
gse125741_gene_name<-gene_ids[["GSE125741"]][!(str_detect(gene_ids[["GSE125741"]],"SAUSA"))]
gse130777_USA<-gene_ids[["GSE130777"]][str_detect(gene_ids[["GSE130777"]],"SAUSA")]
gse130777_gene_name<-gene_ids[["GSE130777"]][!(str_detect(gene_ids[["GSE130777"]],"SAUSA"))]
gse132179_USA<-gene_ids[["GSE132179"]]
gse148024_USA<-unlist(lapply(gene_ids[["GSE148024"]],function(x){str_extract(x,'SAUSA300\\_[0-9]+')}))

#GSE202268 Newman study added 1/17/24
newman_RS_table<-read_xlsx("/Users/lilytaub/Downloads/GeneSpecificInformation_Newman.xlsx", col_names = T)
length(gene_ids[["GSE202268"]])==sum(str_detect(gene_ids[["GSE202268"]], "NWMN"))
gse202268_NWM<-NWMN_RS_dict[NWMN_RS_dict$`new locus tag` %in% gene_ids[["GSE202268"]],]$`locus tag`

#GSE189644 USA300 study added 1/18/24
gse189644_USA<-str_replace(gene_ids[["GSE189644"]], pattern="gene\\.","")
length(gse189644_USA)==sum(str_detect(gse189644_USA, "USA"))
gse189644_USA<-USA300_RS_dict[USA300_RS_dict$`new locus tag` %in% gse189644_USA,]$`locus tag`

#GSE239718 USA300 study added 1/18/24
gse239718_gene_name<-gene_ids[["GSE239718"]][which(!(str_detect(gene_ids[["GSE239718"]],"USA")))]
gse239718_USA<-gene_ids[["GSE239718"]][which(str_detect(gene_ids[["GSE239718"]],"USA"))]
gse239718_USA <- USA300_RS_dict[USA300_RS_dict$`new locus tag` %in% gse239718_USA,]$`locus tag`

USA_Numbers<-union(union(union(union(union(gse125741_USA,gse130777_USA),gse132179_USA),gse148024_USA),gse189644_USA),gse239718_USA)
USA_gene_names<-union(union(gse125741_gene_name,gse130777_gene_name),gse239718_gene_name)
NWMN_Numbers <- union(gse122048_NWM,gse202268_NWM)

write.csv(NWMN_Numbers,file=paste(outdir,"NWM_Numbers.txt",sep=""),row.names=F,quote=F)
write.csv(gse122048_gene_name,file=paste(outdir,"NWM_gene_names.txt",sep=""),row.names = F,quote = F)
write.csv(USA_Numbers,file=paste(outdir,"USA_Numbers.txt",sep=""),row.names = F,quote=F)
write.csv(USA_gene_names,file=paste(outdir,"USA_gene_names.txt",sep=""),row.names = F,quote = F)

# Streptococcus_Species
gse89964_SSA<-gene_ids[["GSE89964"]][str_detect(gene_ids[["GSE89964"]],"SSA")]
gse90021_SSA<-gene_ids[["GSE90021"]][str_detect(gene_ids[["GSE90021"]],"SSA")]
gse90021_gene_name<-gene_ids[["GSE90021"]][!(str_detect(gene_ids[["GSE90021"]],"SSA"))]
gse97218_SSA<-gene_ids[["GSE97218"]][str_detect(gene_ids[["GSE97218"]],"SSA")]
gse97357_SSA<-gene_ids[["GSE97357"]][str_detect(gene_ids[["GSE97357"]],"SSA")]
# GSE209925 added 1/18/24
gse209925_SSA <- SK36_RS_dict[SK36_RS_dict$New_Locus_Tag %in% gene_ids[["GSE209925"]],]$Locus_Tag

SSA_Numbers<-union(union(union(union(gse89964_SSA,gse90021_SSA),gse97218_SSA),gse97357_SSA),gse209925_SSA)

write.csv(SSA_Numbers,file=paste(outdir,"SSA_Numbers.txt",sep=""),row.names = F,quote=F)
write.csv(gse90021_gene_name,file=paste(outdir,"SSA_gene_names.txt",sep=""),row.names = F,quote = F)


# UniProt Translation -----------------------------------------------------

# files with Uniprot translations, translated with https://www.uniprot.org/id-mapping
files<-c("BT_Numbers.tsv","NWM_gene_names.tsv","NWM_Numbers.tsv","PA14_gene_names.tsv","PA14_Numbers.tsv",
         "PAO1_gene_names.tsv","PAO1_Numbers.tsv","SSA_gene_names.tsv","SSA_Numbers.tsv",
         "USA300_gene_names.tsv","USA300_Numbers.tsv") 

Gene_ID_Dictionary<-data.frame("From"=NA,"Entry"=NA)
strain_dict<-list("BT"="VPI-5482","NWM"="Newman","PA14"="PA14","PAO1"="PAO1","SSA"="SK36","USA300"="USA300")

for (f in files){
  translation<-read.delim(file=paste(transdir,f,sep=""))
  translation$Strain<-strain_dict[[unlist(str_split(f,"_"))[1]]]
  Gene_ID_Dictionary<-bind_rows(Gene_ID_Dictionary,translation)
}
Gene_ID_Dictionary<-na.omit(Gene_ID_Dictionary)

#-- add NWM, USA300, and SK36 data from CF-Seq UniProt table to Gene_ID_Dictionary
{
  idx_to_add<-which(!(NWMN_UniProt_dict$GENENAME %in% Gene_ID_Dictionary[Gene_ID_Dictionary$Strain=="Newman",]$From))
  lookup<-c(From="GENENAME",Entry="UNIPROTKB")
  NWMN_to_add<-NWMN_UniProt_dict[idx_to_add,] %>% 
    select(c(UNIPROTKB,GENENAME)) %>% 
    rename(all_of(lookup)) %>% 
    mutate(Strain="Newman")
  Gene_ID_Dictionary <- bind_rows(Gene_ID_Dictionary,NWMN_to_add)
  
  idx_to_add <- which(!(USA300_UniProt_dict$GENENAME %in% Gene_ID_Dictionary[Gene_ID_Dictionary$Strain=="USA300",]$From))
  USA300_to_add <- USA300_UniProt_dict[idx_to_add,] %>% 
    rename(all_of(lookup)) %>% 
    mutate(Strain="USA300")
  Gene_ID_Dictionary <- bind_rows(Gene_ID_Dictionary, USA300_to_add)
  
  idx_to_add <- which(!(SK36_UniProt_dict$GENENAME %in% Gene_ID_Dictionary[Gene_ID_Dictionary$Strain=="SK36",]$From))
  SK36_to_add <- SK36_UniProt_dict[idx_to_add,] %>% 
    rename(all_of(lookup)) %>% 
    mutate(Strain="SK36")
  Gene_ID_Dictionary <- bind_rows(Gene_ID_Dictionary, SK36_to_add)
}

# non-standard format studies:
skip<-c("GSE152480","GSE153067","GSE163234","GSE166986","GSE89791")
# GSE89791 --> Inter pattern (didn't translate to uniprot)

# remove non-standard format studies from Raw_Count_Data
Raw_Count_Data<-Raw_Count_Data[which(!names(Raw_Count_Data) %in% skip)]


# add Uniprot IDS as list column
App_Data_Uniprot<-App_Data

#-- Dropped Studies --#
{
  #GSE153067 --> not unique study comparisons, dropped
  for (comp in names(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE153067"]][-(1)])){
    print(comp)
    App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE153067"]][[comp]]<-App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE153067"]][[comp]][-(3251),]
    rownames(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE153067"]][[comp]])<-gse153067_PA[-(3251)]
  }

  #GSE166986 --> not unique rownames, dropped
  from<-data.frame("from"=rownames(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE166986"]][["Delta.rhlE1 - Delta.rhlE1...Delta.rhlE2 Full Results"]]))
  from<-from %>% 
    left_join(gse166986,by=join_by("from"=="CDS number"))
  for (comp in names(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE166986"]][-(1)])){
    rownames(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE166986"]][[comp]])<-from$`PA number` 
  }

  # GSE152480 --> non unique rownames when extracting PA number, dropped
  pa_num_idx<-which(str_detect(gene_ids[["GSE152480"]],"PA")==TRUE)
  for (comp in names(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE152480"]][-(1)])){
    rownames(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE152480"]][[comp]][pa_num_idx,])<-gse152480_PA[pa_num_idx]
  }
}

#-- Special Cases --#
{
  # GSE148024
  to<-str_replace(rownames(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE148024"]][[4]]),"gene\\.","")
  App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE148024"]][["Additional_Metadata"]]$Strain<-"USA300, Cowan, Clinical isolates"
  strains[["GSE148024"]]<-"USA300, Cowan, Clinical isolates"
  for (comp in names(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE148024"]][-(1)])){
    rownames(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE148024"]][[comp]])<-to
  }
  # GSE124385
  from<-data.frame("from"=rownames(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE124385"]][["X.rclR.Control. - X.rclR.HOCl Full Results"]]))
  RS_idx<-which(str_detect(from$from,"RS"))
  from<-from %>% 
    left_join(RS_gene_dict[str_detect(RS_gene_dict$new_locus_tag,"PA14"),],by=join_by("from"=="old_locus_tag"))
  to<-from$new_locus_tag
  NA_idx<-which(is.na(to==TRUE))
  to[NA_idx]<-from$from[NA_idx]
  for (comp in names(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE124385"]][-(1)])){
    rownames(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE124385"]][[comp]])<-to
  }
  # GSE142464
  from<-data.frame("from"=rownames(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE142464"]][["B136.33.1 - C3719.3 Full Results"]]))
  RS_idx<-which(str_detect(from$from,"RS"))
  from<-from %>% 
    left_join(RS_gene_dict[str_detect(RS_gene_dict$new_locus_tag,"PA14"),],by=join_by("from"=="old_locus_tag"))
  to<-from$new_locus_tag
  NA_idx<-which(is.na(to==TRUE))
  to[NA_idx]<-from$from[NA_idx]
  for (comp in names(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE142464"]][-(1)])){
    rownames(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE142464"]][[comp]])<-to
  }
  #GSE185398
  from<-data.frame("from"=rownames(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE185398"]][[2]]))
  RS_idx<-which(str_detect(from$from,"RS"))
  from<-from %>% 
    left_join(RS_gene_dict[str_detect(RS_gene_dict$new_locus_tag,"PA14"),],by=join_by("from"=="old_locus_tag"))
  to<-from$new_locus_tag
  NA_idx<-which(is.na(to==TRUE))
  to[NA_idx]<-from$from[NA_idx]
  for (comp in names(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE185398"]][-(1)])){
    rownames(App_Data_Uniprot[["Pseudomonas_aeruginosa"]][["GSE185398"]][[comp]])<-to
  }
  #GSE202268
  from<-data.frame("from"=rownames(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE202268"]][[2]]))
  RS_idx<-which(str_detect(from$from,"RS"))
  from<-from %>% 
    left_join(NWMN_RS_dict[str_detect(NWMN_RS_dict$`new locus tag`,"NWMN"),],by=join_by("from"=="new locus tag")) %>% 
    select(-UniProt)
  to<-from$`locus tag`
  NA_idx<-which(is.na(to==TRUE))
  to[NA_idx]<-from$from[NA_idx]
  for (comp in names(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE202268"]][-(1)])){
    rownames(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE202268"]][[comp]])<-to
  }
  #GSE189644
  from<-data.frame("from"=rownames(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE189644"]][[2]])) %>% 
    mutate(from=str_replace(from, "gene\\.",""))
  RS_idx<-which(str_detect(from$from,"RS"))
  from<-from %>% 
    left_join(USA300_RS_dict[str_detect(USA300_RS_dict$`new locus tag`,"USA300"),],by=join_by("from"=="new locus tag")) %>% 
    select(-UniProt)
  to<-from$`locus tag`
  NA_idx<-which(is.na(to==TRUE))
  to[NA_idx]<-from$from[NA_idx]
  for (comp in names(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE189644"]][-(1)])){
    rownames(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE189644"]][[comp]])<-to
  }
  #GSE239718
  from<-data.frame("from"=rownames(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE239718"]][[2]]))
  RS_idx<-which(str_detect(from$from,"RS"))
  from<-from %>% 
    left_join(USA300_RS_dict[str_detect(USA300_RS_dict$`new locus tag`,"USA300"),],by=join_by("from"=="new locus tag")) %>% 
    select(-UniProt)
  to<-from$`locus tag`
  NA_idx<-which(is.na(to==TRUE))
  to[NA_idx]<-from$from[NA_idx]
  for (comp in names(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE239718"]][-(1)])){
    rownames(App_Data_Uniprot[["Staphylococcus_aureus"]][["GSE239718"]][[comp]])<-to
  }
  #GSE209925
  from<-data.frame("from"=rownames(App_Data_Uniprot[["Streptococcus_Species"]][["GSE209925"]][[2]]))
  RS_idx<-which(str_detect(from$from,"RS"))
  from<-from %>% 
    left_join(SK36_RS_dict,by=join_by("from"=="New_Locus_Tag"))
  to<-from$Locus_Tag
  NA_idx<-which(is.na(to==TRUE))
  to[NA_idx]<-from$from[NA_idx]
  for (comp in names(App_Data_Uniprot[["Streptococcus_Species"]][["GSE209925"]][-(1)])){
    rownames(App_Data_Uniprot[["Streptococcus_Species"]][["GSE209925"]][[comp]])<-to
  }
}

# translate genes to Uniprot
for (species in names(App_Data_Uniprot)){
  for (study in names(App_Data_Uniprot[[species]])){
    if (study %in% skip){
      next
    }
    comp<-names(App_Data_Uniprot[[species]][[study]])[2]
    genes<-rownames(App_Data_Uniprot[[species]][[study]][[comp]])
    strain<-unlist(str_split(strains[[study]],","))[1]
    uniprot<-list()
    print(species)
    print(study)
    for (i in seq(1,length(genes))){
      if (genes[i] %in% Gene_ID_Dictionary[Gene_ID_Dictionary$Strain==strain,]$From){
        to_add<-as.list(Gene_ID_Dictionary[Gene_ID_Dictionary$From==genes[i] & Gene_ID_Dictionary$Strain==strain,]$Entry)
        uniprot[[i]]<-to_add
      }else{
        uniprot[[i]]<-NA
      }
    }
    for (comp in names(App_Data_Uniprot[[species]][[study]])[-(1)]){
      App_Data_Uniprot[[species]][[study]][[comp]]$Uniprot<-uniprot
    }
  }
}

# Get expression data (logFC) per study comparison
eskape_dir<-paste0(outdir,"For_ESKAPE")
For_ESKAPE<-list()

for (species in species_list){
  for (study in names(App_Data_Uniprot[[species]])){
    if (study %in% skip){
      next
    }
    for (comp in names(App_Data_Uniprot[[species]][[study]][-(1)])){
      
      expression_data<-App_Data_Uniprot[[species]][[study]][[comp]] %>% 
        unnest(cols="Uniprot") %>% 
        filter(!(Uniprot=="NULL")) %>% 
        select(c("Uniprot","logFC","PValue"))
      
      expression_data$Uniprot<-unlist(expression_data$Uniprot) 
      filename=paste(study,comp,".csv",sep="_")
      write.csv(expression_data,file=paste(eskape_dir,filename,sep="/"),quote = F)
      For_ESKAPE[[species]][[study]][[comp]]<-expression_data
      For_ESKAPE[[species]][[study]][["Additional_Metadata"]]<-App_Data_Uniprot[[species]][[study]][["Additional_Metadata"]]
    }
  }
}

save(For_ESKAPE,Raw_Count_Data,file=paste(outdir,"For_Eskape.Rdata",sep=""))
