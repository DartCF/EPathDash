## PATHWAY ANALYSIS FOR EPathDash ##

library(stringr)

# directory definitions
datadir<-"/Users/lilytaub/Documents/ESKAPE_Act_Plus_V3/"
datadir2<-"/Users/lilytaub/Documents/DartCF/SummerCapstone/DataCollection/"
keggdir<-"/Users/lilytaub/Documents/DartCF/SummerCapstone/DataCollection/KEGG_Results/"
godir<-"/Users/lilytaub/Documents/DartCF/SummerCapstone/DataCollection/GO_Results/"
appdir<-"/Users/lilytaub/Documents/DartCF/SummerCapstone/PathwayApplication/"

# load in ESKAPE database
load(paste(datadir,"KEGG_GO.Rdata",sep=""))
# load in Study data formatted for ESKAPE analysis (done with data_collection.R script)
load(paste(datadir2,"For_Eskape.Rdata",sep=""))

# subset relevant strains
KEGG<-KEGG[c("pau","pae","bth","saa","sae","ssa")]
GO<-GO[c("pau","pae","bth","saa","sae","ssa")]
Path_List<-Path_List[c("pau","pae","bth","saa","sae","ssa")]


KEGGActivationAnalysis<-function(strain,Results,studyID,comp){
  strain_dict<-list("PA14"="pau","PAO1"="pae","VPI-5482"="bth","USA300"="saa","Newman"="sae","SK36"="ssa")
  strain<-strain_dict[[strain]]
  print(paste("KEGG",strain,sep=":"))
  KEGG_Prots <- KEGG[[strain]]$Entry[KEGG[[strain]]$Entry %in% Results$V1]
  
  # select paths from KEGG strain that match genes in input file
  KEGG_Paths <- KEGG[[strain]]$KEGG_paths[KEGG[[strain]]$Entry %in% Results$V1]
  
  # turning table back into list:
  # list of genes involved in each path
  KEGG_List <- tapply(KEGG_Prots, KEGG_Paths, c)
  
  KEGG_List_Length <- lapply(KEGG_List, length)
  Good_KEGG <- names(KEGG_List)[KEGG_List_Length > 4]
  
  KEGG_Tests <- lapply(Good_KEGG, function (kname){ #kname is pathway name
    x = KEGG_List[[kname]] # genes involved in path
    n = sum(unique(KEGG_Prots) %in% x) # number of KEGG genes in genes in KEGG path
    if (n > 0) {
      #gene expression values for genes from user input data which are in path
      o = Results[which(unique(Results$V1) %in% x), "V2"] 
      s = sum (o > 0) #number of genes with > 0 expression
      #likelihood that probability of expression>0 given number of positively expressed genes and the number of genes in the path
      binom.test(x = s, p = 0.5, n = n) 
    } 
    else {
      binom.test(x = 1, p = 0.5, n = 1)
    }
    
  })
  names(KEGG_Tests) <- Good_KEGG
  
  P <- unlist(lapply(KEGG_Tests, function(x){x$p.value}))
  Est <- unlist(lapply(KEGG_Tests, function(x){x$estimate}))
  FDR <- p.adjust(P, "fdr") #adjusting for testing multiple hypothesis (Benjamini & Hochberg)
  
  KEGG_FC <- lapply(KEGG_List[Good_KEGG], function(x){
    Results[unique(Results$V1) %in% x, "V2"]
  })
  
  KEGG_FC_median <- lapply(KEGG_FC, median)
  
  KEGG_Result <- data.frame(
    "Path" = Good_KEGG,
    "Estimate" = round(as.numeric(Est), 2),
    "Median_fold_change" = round(as.numeric(KEGG_FC_median), 2),
    "P_value" = as.numeric(P),
    "FDR" = as.numeric(FDR))
  filename<-paste(studyID,comp,sep="_")
  write.csv(KEGG_Result, paste(keggdir,filename,".csv",sep=""), row.names = FALSE)
  
  Sig <- KEGG_Result
  Sig <- Sig[Sig$FDR < 0.05, ]
  
  KEGGSigList <- KEGG_List[Sig$Path]

  return(KEGGSigList)
}

GOActivationAnalysis<-function(strain,Results,studyID,comp){
  
  strain_dict<-list("PA14"="pau","PAO1"="pae","VPI-5482"="bth","USA300"="saa","Newman"="sae","SK36"="ssa")
  strain<-strain_dict[[strain]]
  
  GO_Prots <- GO[[strain]]$Entry[GO[[strain]]$Entry %in% Results$V1]
  GO_Paths <- GO[[strain]]$Gene.ontology..GO.[GO[[strain]]$Entry %in% Results$V1]
  
  # turning table back into list:
  GO_List <- tapply(GO_Prots, GO_Paths, c)
  
  # require at least 4 genes per path:
  GO_List_Length <- lapply(GO_List, length)
  Good_GO <- names(GO_List)[GO_List_Length > 4]
  
  GO_Tests <- lapply(Good_GO, function (kname){
    x = GO_List[[kname]]
    n = sum(unique(GO_Prots) %in% x)
    if (n > 0) {
      o = Results[which(unique(Results$V1) %in% x), "V2"]
      s = sum (o > 0)
      binom.test(x = s, p = 0.5, n = n)    
    } else {binom.test(x = 1, p = 0.5, n = 1)}
  })
  names(GO_Tests) <- Good_GO
  
  GO_Tests_P <- unlist(lapply(GO_Tests, function(x){x$p.value}))
  
  GO_Tests_Est <- unlist(lapply(GO_Tests, function(x){x$estimate}))
  names(GO_Tests_Est) <- names(GO_Tests_P)
  
  GO_FDR <- p.adjust(GO_Tests_P, "fdr")
  
  GO_FC <- lapply(GO_List[Good_GO], function(x){
    Results[which(unique(Results$V1) %in% x), "V2"]
  })
  
  GO_FC_median <- lapply(GO_FC, median)
  
  
  GO_Result <- data.frame(
    "GO_term" = names(GO_Tests),
    "Estimate" = round(as.numeric(GO_Tests_Est), 2),
    "Median_fold_change" = round(as.numeric(GO_FC_median), 2),
    "P_value" = as.numeric(GO_Tests_P),
    "FDR" = as.numeric(GO_FDR))
  filename<-paste(studyID,comp,sep="_")
  write.csv(GO_Result, paste(godir,filename,".csv",sep=""), row.names = FALSE)
  
  
  GO_Sig <- GO_Result[GO_Result$FDR < 0.05, ]
  GO_SigList <- GO_List[GO_Sig$GO_term]
  return(GO_SigList)
}


# Run Pathway Activation Analysis -----------------------------------------

# Sig_Lists: for each comparison within a study, list activated/repressed pathways for KEGG and GO
ESKAPE_Sig_Lists<-list()
for (species in names(For_ESKAPE)){
  for (study in names(For_ESKAPE[[species]])){
    strain<-unlist(str_split(For_ESKAPE[[species]][[study]][["Additional_Metadata"]][["Strain"]],","))[1]
    for (comp in names(For_ESKAPE[[species]][[study]])){
      if (str_detect(comp,"Additional")){
        next
      }else{
        input_data<-data.frame(For_ESKAPE[[species]][[study]][[comp]])
        colnames(input_data)<-c("V1","V2")
        print(study)
        comp<-str_split_1(comp," Full Results")[1]
        ESKAPE_Sig_Lists[[species]][[study]][[comp]][["KEGG"]]<-KEGGActivationAnalysis(strain,input_data,study,comp)
        ESKAPE_Sig_Lists[[species]][[study]][[comp]][["GO"]]<-GOActivationAnalysis(strain,input_data,study,comp)
      }
    }
  }
}

# Read results into a data structure (all results, not just significant)
ESKAPE_Results<-list()
for (species in names(ESKAPE_Sig_Lists)){
  for (study in names(ESKAPE_Sig_Lists[[species]])){
    for (comp in names(ESKAPE_Sig_Lists[[species]][[study]])){
      filename<-paste(study,comp,sep="_")
      results_kegg<-read.csv(paste(keggdir,filename,".csv",sep=""))
      results_go<-read.csv(paste(godir,filename,".csv",sep=""))
      ESKAPE_Results[[species]][[study]][[comp]][["KEGG"]]<-results_kegg
      ESKAPE_Results[[species]][[study]][[comp]][["GO"]]<-results_go
    }
  }
}

# create metadata only data object for app
Study_Metadata<-list()
Study_LogFC_Genes<-list()
for (species in names(For_ESKAPE)){
  for (study in names(For_ESKAPE[[species]])){
    Study_Metadata[[species]][[study]]<-For_ESKAPE[[species]][[study]][["Additional_Metadata"]]
    Study_Metadata[[species]][[study]]["LinkNCBI"]<-paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc",study,sep="=")
    new_colnames<-append(c("Date","GEO.Accession","Title"),colnames(Study_Metadata[[species]][[study]])[4:ncol(Study_Metadata[[species]][[study]])])
    colnames(Study_Metadata[[species]][[study]])<-new_colnames
    for(comp in names(For_ESKAPE[[species]][[study]])){
      if (str_detect(comp, "Additional")){
        next
      }else{
        comp_rename<-str_split_1(comp," Full Results")[1]
        Study_LogFC_Genes[[study]][[comp_rename]]<-For_ESKAPE[[species]][[study]][[comp]]
      }
    }
  }
}

# KEGG pathway data for KEGG mapper URL
names(KEGG)<-c("PA14","PAO1","VPI-5482","USA300","Newman","SK36")
names(GO)<-c("PA14","PAO1","VPI-5482","USA300","Newman","SK36")
names(Path_List)<-c("PA14","PAO1","VPI-5482","USA300","Newman","SK36")

# Rdata file for EPathDash
save(ESKAPE_Results,ESKAPE_Sig_Lists,Study_Metadata,Study_LogFC_Genes,KEGG,GO,Path_List,Raw_Count_Data,file=paste(appdir,"Data.Rdata",sep=""))
