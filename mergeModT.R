

# This function is designed to parse selected columns of individual ModT table outputs and merge by using Uniprot IDs

# By default this function filters out keratins from the lists


mergeModT<-function(ModT_ClassV_match ) {
  
  
  mt<-read.csv(ModT_ClassV_match,header=TRUE) # Read the ModT_ClassV match table
  
  colnames(mt)<-c("ModT","ClassV","Experiment", "Cell")
  
  n<-nrow(mt)
  
            #####################################################################################################################################
            MT<-read.csv(file=as.character(mt$ModT[1]),header=TRUE,check.names = F) #Read the ModT results table for the first experiment
            ClV<-read.csv(file=as.character(mt$ClassV[1]),header=TRUE,check.names = F) #Read the Class Vector table for the experiment
            Cell<-as.character(mt$Cell[1]) # Name of the cell line
            Expr<-as.character(mt$Experiment[1]) # Name of the Experiment
            ######################################################################################################################################
            w<-which(MT$species!="HUMAN") #Looking for any non-human contaminants
            if(length(w)>0) {MT<-MT[-w,]} #Removes any non-human contaminants
            
            KRT<-read.csv("human_keratins_as_downloaded_from_HGNC_06242016.csv",header =T)
            KRT<-KRT$Approved
            w<-which(MT$geneSymbol %in%  KRT)
            if(length(w)>0) {MT<-MT[-w,]} #Removes any keratins            
            #######################################################################################################################################
            
            ##Exract treatment-specific information from this ModT table
            
                    a<-1
                    n1<-nrow(ClV)/2
                    
                    Edf<-data.frame(MT[,"accession_number"],check.names = F)
                    colnames(Edf)<-c("accession_number")
                    
                    for(i in 1:n1 ) {
                      
                      
                      
                      adjp<-which(colnames(MT)==paste("adj.P.Val",ClV[a,2],sep=".")) #Extract relevant column indexes for the specific treatment
                      AExp<-which(colnames(MT)==paste("AveExpr",ClV[a,2],sep="."))
                      p<-which(colnames(MT)==paste("P.Value",ClV[a,2],sep="."))
                      
                      
                      Etemp<-data.frame(MT[,AExp],MT[,adjp],MT[,p])#Extract the relevant columns for a specified treatment and coerce them into a data frame
                      
                          colnames(Etemp)<-c(paste("Average_Fold_Change(log2)",Cell,ClV[a,2],sep="_"),
                                           paste("adj.P.Val",Cell,ClV[a,2],sep="_"),
                                           paste("P.Value",Cell,ClV[a,2],sep="_"))
                  
                          Edf<-data.frame(Edf,Etemp,check.names = F)    
                          
                          a<-a+2 #Move to the next treatment within the experiment, until no treatment left
                          
                    } 
                    
                    
                    #Adding more general columns for this ModT table
                    
                    ##geneSymbol column
                    
                    w<-which(colnames(MT) == "geneSymbol")
                    
                    Edf<-data.frame(Edf, MT[,w], check.names = F)
                    colnames(Edf)[ncol(Edf)]<-paste("geneSymbol",Cell,Expr,sep = "_")
                    
                    
                    ##TMT ratio columns
                    
                    w<-which(colnames(MT) %in%  ClV[,1]) # Find the columns in the MT that represent TMT ratios
                    
                    Edf<-data.frame(Edf, MT[,w], check.names = F)
                    
                    ##Unique peptides 
                    
                    NAME<-sub('\\..*', '', ClV[1,1]) #Get the directory name for the experiment from the class vector (anything before the first dot ".")
                    
                    w<-which(colnames(MT)==paste(NAME,"unique_peptides",sep="."))
                    Edf<-data.frame(Edf, MT[,w])
                    colnames(Edf)[ncol(Edf)]<-colnames(MT)[w]
                    
                    
                    ##Total intensity
                    w<-which(colnames(MT)==paste(NAME,"totalIntensity",sep="."))
                    Edf<-data.frame(Edf, MT[,w])
                    colnames(Edf)[ncol(Edf)]<-colnames(MT)[w]
                    
                    
                                    for (i in 2:n) { #Loop over the rest of the experiments specified in the match table
                          
                                      #####################################################################################################################################
                                      MT<-read.csv(file=as.character(mt$ModT[i]),header=TRUE,check.names = F) #Read the ModT results table for the first experiment
                                      ClV<-read.csv(file=as.character(mt$ClassV[i]),header=TRUE,check.names = F) #Read the Class Vector table for the experiment
                                      Cell<-as.character(mt$Cell[i]) # Name of the cell line
                                      Expr<-as.character(mt$Experiment[i]) # Name of the Experiment
                                      ######################################################################################################################################
                                      w<-which(MT$species!="HUMAN") #Looking for any non-human contaminants
                                      if(length(w)>0) {MT<-MT[-w,]} #Removes any non-human contaminants
                                      
                                      KRT<-read.csv("human_keratins_as_downloaded_from_HGNC_06242016.csv",header =T)
                                      KRT<-KRT$Approved
                                      w<-which(MT$geneSymbol %in%  KRT)
                                      if(length(w)>0) {MT<-MT[-w,]} #Removes any keratins            
                                      #######################################################################################################################################
                                      
                                      ##Exract treatment-specific information from this ModT table
                                      
                                      a<-1
                                      n1<-nrow(ClV)/2
                                      
                                      Edftemp<-data.frame(MT[,"accession_number"],check.names = F)
                                      colnames(Edftemp)<-c("accession_number")
                                      
                                                    for(i in 1:n1 ) {
                                                      
                                                      
                                                      
                                                      adjp<-which(colnames(MT)==paste("adj.P.Val",ClV[a,2],sep=".")) #Extract relevant column indexes for the specific treatment
                                                      AExp<-which(colnames(MT)==paste("AveExpr",ClV[a,2],sep="."))
                                                      p<-which(colnames(MT)==paste("P.Value",ClV[a,2],sep="."))
                                                      
                                                      
                                                      Etemp<-data.frame(MT[,AExp],MT[,adjp],MT[,p])#Extract the relevant columns for a specified treatment and coerce them into a data frame
                                                      
                                                      colnames(Etemp)<-c(paste("Average_Fold_Change(log2)",Cell,ClV[a,2],sep="_"),
                                                                         paste("adj.P.Val",Cell,ClV[a,2],sep="_"),
                                                                         paste("P.Value",Cell,ClV[a,2],sep="_"))
                                                      
                                                      Edftemp<-data.frame(Edftemp,Etemp,check.names = F)    
                                                      
                                                      a<-a+2 #Move to the next treatment within the experiment, until no treatment left
                                                      
                                                    } 
                                      
                                      
                                                      #Adding more general columns for this ModT table
                                                      
                                                      ##geneSymbol column
                                                      
                                                      w<-which(colnames(MT) == "geneSymbol")
                                                      
                                                      Edftemp<-data.frame(Edftemp, MT[,w], check.names = F)
                                                      colnames(Edftemp)[ncol(Edftemp)]<-paste("geneSymbol",Cell,Expr,sep = "_")
                                                      
                                                      
                                                      ##TMT ratio columns
                                                      
                                                      w<-which(colnames(MT) %in%  ClV[,1]) # Find the columns in the MT that represent TMT ratios
                                                      
                                                      Edftemp<-data.frame(Edftemp, MT[,w], check.names = F)
                                                      
                                                      ##Unique peptides 
                                                      
                                                      NAME<-sub('\\..*', '', ClV[1,1]) #Get the directory name for the experiment from the class vector (anything before the first dot ".")
                                                      
                                                      w<-which(colnames(MT)==paste(NAME,"unique_peptides",sep="."))
                                                      Edftemp<-data.frame(Edftemp, MT[,w])
                                                      colnames(Edftemp)[ncol(Edftemp)]<-colnames(MT)[w]
                                                      
                                                      
                                                      ##Total intensity
                                                      w<-which(colnames(MT)==paste(NAME,"totalIntensity",sep="."))
                                                      Edftemp<-data.frame(Edftemp, MT[,w])
                                                      colnames(Edftemp)[ncol(Edftemp)]<-colnames(MT)[w]
                                      
                                      
                                      
                                      ##Merge Individual ModTtables here
                                      
                                      
                                      Edf<-merge(Edf, Edftemp, by.x="accession_number", by.y ="accession_number", all = T )
                          
                                    
                                       # Continue until all the experiments listed in the match table are completed                                    
                                      } 
  
                    write.csv(file="Merged_datasets.csv",Edf)
  
}