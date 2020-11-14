RNEA<-function(filename,identifier="GeneName",species,FC_threshold=1,PV_threshold=0.05,output="Output",network="regulatory",type_of_output="html"){
  options(warn=-1);
  if(!require(SortableHTMLTables)){
    install.packages("SortableHTMLTables");
  }
  library(SortableHTMLTables);
  if(identifier=="Refseq"){
    Refseq=read.table("RNEA/RNEA/ReferenceFiles/Refseq2Gene.txt",row.names=2,sep="\t");
  }
  TF_ref=paste("RNEA/RNEA/ReferenceFiles/",species,"_TF_plus_oreganno_trrust.tsv",sep="");
  
  kegg_ref=paste("RNEA/RNEA/ReferenceFiles/",species,"_kegg.tsv",sep="");
  keggcat_ref=paste("RNEA/RNEA/ReferenceFiles/",species,"_keggcat.tsv",sep="");
  GOs_ref=paste("RNEA/RNEA/ReferenceFiles/",species,"_GOs.tsv",sep="");
  
  TF=read.table(TF_ref,sep="\t",fill=T);
  rownames(TF)=TF[,2];	
  
  kegg=read.table(kegg_ref,sep="\t",fill=T);
  rownames(kegg)=kegg[,2];
  keggcat=read.table(keggcat_ref,sep="\t",fill=T);
  rownames(keggcat)=keggcat[,2];
  gos=read.table(GOs_ref,sep="\t",fill=T);
  rownames(gos)=gos[,2];
  Dereg_Targets=list();
  
  Input=read.table(filename,sep="\t",header=T);
  if(ncol(Input)==3){
    Total_number_of_genes=length(unique(Input[,1]));
  }else if(ncol(Input)==14){
    Total_number_of_genes=length(unique(Input[,3]));
  }
  TF_counts=cbind(TF[,1],Total_number_of_genes,0,0,0,0);
  colnames(TF_counts)=c("Targets","Genes","IsDE?","TargetsDE","TargetsUp","TargetsDown");
  rownames(TF_counts)=TF[,2];
  
  kegg_counts=cbind(kegg[,1],Total_number_of_genes,0,0,0);
  colnames(kegg_counts)=c("Members","Genes","MembersDE","MembersUp","MembersDown");
  rownames(kegg_counts)=kegg[,2];
  keggcat_counts=cbind(keggcat[,1],Total_number_of_genes,0,0,0);
  colnames(keggcat_counts)=c("Members","Genes","MembersDE","MembersUp","MembersDown");
  rownames(keggcat_counts)=keggcat[,2];
  go_counts=cbind(gos[,1],Total_number_of_genes,0,0,0);
  colnames(go_counts)=c("Members","Genes","MembersDE","MembersUp","MembersDown");
  rownames(go_counts)=gos[,2];
  Number_of_DE_genes=0;	
  Number_of_Up_genes=0;
  Number_of_Down_genes=0;
  DE_genes=matrix(data=0,nrow=1);
  UP_genes=matrix(data=0,nrow=1);
  DOWN_genes=matrix(data=0,nrow=1);
  Network=matrix(ncol=2,data="");
  colnames(Network)=c("Source","Target");
  
  
  if(ncol(Input)==14){
    for(line in 2:nrow(Input)){
      if((Input[line,13]<=PV_threshold)&&(Input[line,10]>=FC_threshold)){
        Number_of_DE_genes=Number_of_DE_genes+1;
        Number_of_Up_genes=Number_of_Up_genes+1;
        if(identifier=="GeneName"){
          DE_genes[Number_of_DE_genes]=as.character(Input[line,3]);
          UP_genes[Number_of_Up_genes]=as.character(Input[line,3]);
        }else if(identifier=="Refseq"){
          DE_genes[Number_of_DE_genes]=as.character(Refseq[Input[line,3],1]);
          UP_genes[Number_of_Up_genes]=as.character(Refseq[Input[line,3],1]);
        }else{
          print("Please select a supported gene identifier (\"GeneName\" or \"Refseq\")");
        }
      }
      if((Input[line,13]<=PV_threshold)&&(Input[line,10]<=-FC_threshold)){
        Number_of_DE_genes=Number_of_DE_genes+1;
        Number_of_Down_genes=Number_of_Down_genes+1;
        if(identifier=="GeneName"){
          DE_genes[Number_of_DE_genes]=as.character(Input[line,3]);
          DOWN_genes[Number_of_Down_genes]=as.character(Input[line,3]);
        }else if(identifier=="Refseq"){
          DE_genes[Number_of_DE_genes]=as.character(Refseq[Input[line,3],1]);
          DOWN_genes[Number_of_Down_genes]=as.character(Refseq[Input[line,3],1]);
        }
        else{
          print("Please select a supported gene identifier (\"GeneName\" or \"Refseq\")");
        }
      }
    }
  }
  if(ncol(Input)==3){
    for(line in 2:nrow(Input)){
      if((Input[line,3]<=PV_threshold)&&(Input[line,2]>=FC_threshold)){
        Number_of_DE_genes=Number_of_DE_genes+1;
        Number_of_Up_genes=Number_of_Up_genes+1;
        if(identifier=="GeneName"){
          DE_genes[Number_of_DE_genes]=as.character(Input[line,1]);
          UP_genes[Number_of_Up_genes]=as.character(Input[line,1]);
        }else if(identifier=="Refseq"){
          DE_genes[Number_of_DE_genes]=as.character(Refseq[Input[line,1],1]);
          UP_genes[Number_of_Up_genes]=as.character(Refseq[Input[line,1],1]);
        }else{
          print("Please select a supported gene identifier (\"GeneName\" or \"Refseq\")");
        }
      }
      if((Input[line,3]<=PV_threshold)&&(Input[line,2]<=-FC_threshold)){
        Number_of_DE_genes=Number_of_DE_genes+1;
        Number_of_Down_genes=Number_of_Down_genes+1;
        if(identifier=="GeneName"){
          DE_genes[Number_of_DE_genes]=as.character(Input[line,1]);
          DOWN_genes[Number_of_Down_genes]=as.character(Input[line,1]);
        }else if(identifier=="Refseq"){
          DE_genes[Number_of_DE_genes]=as.character(Refseq[Input[line,1],1]);
          DOWN_genes[Number_of_Down_genes]=as.character(Refseq[Input[line,1],1]);
        }
        else{
          print("Please select a supported gene identifier (\"GeneName\" or \"Refseq\")");
        }
      }
    }
  }
  
  DE_genes=unique(DE_genes);
  UP_genes=unique(UP_genes);
  DOWN_genes=unique(DOWN_genes);
  Number_of_DE_genes=length(DE_genes);
  Number_of_Up_genes=length(UP_genes);
  Number_of_Down_genes=length(DOWN_genes);
  print(paste(Number_of_DE_genes,"differentially expressed genes"));
  print(paste(Number_of_Up_genes,"upregulated"));
  print(paste(Number_of_Down_genes,"downregulated"));
  
  for(Tf in 1:nrow(TF)){
    z<-intersect(TF[Tf,2],UP_genes);
    if(length(z)>0){
      deup<-intersect(as.matrix(TF[z,3:(2+TF[z,1])]),as.list(UP_genes));
      dedown<-intersect(as.matrix(TF[z,3:(2+TF[z,1])]),as.list(DOWN_genes));
      if((network=="global")||(network=="regulatory")){
        if(length(deup>0)){
          for(f in 1:length(deup)){
            Network=rbind(Network,c(z,deup[f]));
          }
        }
        if(length(dedown>0)){
          for(f in 1:length(dedown)){
            Network=rbind(Network,c(z,dedown[f]));
          }
        }
      }
      TF_counts[z,3]=1;
      TF_counts[z,4]=TF_counts[z,4]+length(deup)+length(dedown);
      TF_counts[z,5]=TF_counts[z,5]+length(deup);
      TF_counts[z,6]=TF_counts[z,6]+length(dedown);
      for(tar in 3:(TF[z,1]+2)){
        if(length(intersect(TF[z,tar],rownames(TF)))>0){
          deup<-intersect(as.matrix(TF[as.character(TF[z,tar]),3:(2+TF[as.character(TF[z,tar]),1])]),UP_genes);
          dedown<-intersect(as.matrix(TF[as.character(TF[z,tar]),3:(2+TF[as.character(TF[z,tar]),1])]),DOWN_genes);
          if(length(deup>0)){
            for(f in 1:length(deup)){
              Network=rbind(Network,c(z,as.character(TF[z,tar])));
              Network=rbind(Network,c(as.character(TF[z,tar]),deup[f]));
            }
          }
          if(length(dedown>0)){
            for(f in 1:length(dedown)){
              Network=rbind(Network,c(z,as.character(TF[z,tar])));
              Network=rbind(Network,c(as.character(TF[z,tar]),dedown[f]));
            }
          }
          TF_counts[z,1]=TF_counts[z,1]+TF[as.character(TF[z,tar]),1];
          TF_counts[z,4]=TF_counts[z,4]+length(deup)+length(dedown);
          TF_counts[z,5]=TF_counts[z,5]+length(deup);
          TF_counts[z,6]=TF_counts[z,6]+length(dedown);
          
        }
      }
      
    }
  }
  
  
  for(Tf in 1:nrow(TF)){
    z<-intersect(TF[Tf,2],DOWN_genes);
    if(length(z)>0){
      dedown<-intersect(as.matrix(TF[z,3:(2+TF[z,1])]),as.list(DOWN_genes));
      deup<-intersect(as.matrix(TF[z,3:(2+TF[z,1])]),as.list(UP_genes));
      if((network=="global")||(network=="regulatory")){
        if(length(deup>0)){
          for(f in 1:length(deup)){
            Network=rbind(Network,c(z,deup[f]));
          }
        }
        if(length(dedown>0)){
          for(f in 1:length(dedown)){
            Network=rbind(Network,c(z,dedown[f]));
          }
        }
      }
      TF_counts[z,3]=-1;
      TF_counts[z,4]=TF_counts[z,4]+length(deup)+length(dedown);
      TF_counts[z,6]=TF_counts[z,6]+length(dedown);
      TF_counts[z,5]=TF_counts[z,5]+length(deup);
      for(tar in 3:(TF[z,1]+2)){
        if(length(intersect(TF[z,tar],rownames(TF)))>0){
          deup<-intersect(as.matrix(TF[as.character(TF[z,tar]),3:(2+TF[as.character(TF[z,tar]),1])]),UP_genes);
          dedown<-intersect(as.matrix(TF[as.character(TF[z,tar]),3:(2+TF[as.character(TF[z,tar]),1])]),DOWN_genes);
          if(length(deup>0)){
            for(f in 1:length(deup)){
              Network=rbind(Network,c(z,as.character(TF[z,tar])));
              Network=rbind(Network,c(as.character(TF[z,tar]),deup[f]));
            }
          }
          if(length(dedown>0)){
            for(f in 1:length(dedown)){
              Network=rbind(Network,c(z,as.character(TF[z,tar])));
              Network=rbind(Network,c(as.character(TF[z,tar]),dedown[f]));
            }
          }
          TF_counts[z,1]=TF_counts[z,1]+TF[as.character(TF[z,tar]),1];
          TF_counts[z,4]=TF_counts[z,4]+length(deup)+length(dedown);
          TF_counts[z,5]=TF_counts[z,5]+length(deup);
          TF_counts[z,6]=TF_counts[z,6]+length(dedown);
          
        }
      }
    }
  }
  
  
  
  
  
  
  for(keg in 1:nrow(kegg)){
    z=rownames(kegg)[keg];
    deup<-intersect(as.matrix(kegg[z,3:(2+kegg[z,1])]),as.list(UP_genes));
    dedown<-intersect(as.matrix(kegg[z,3:(2+kegg[z,1])]),as.list(DOWN_genes));
    de<-intersect(as.matrix(kegg[z,3:(2+kegg[z,1])]),as.list(DE_genes));
    if((network=="global")||(network=="functional")){
      if(length(de>0)){
        for(f in 1:length(de)){
          Network=rbind(Network,c(z,de[f]));
        }
      }
    }
    kegg_counts[z,3]=kegg_counts[z,3]+length(de);
    kegg_counts[z,4]=kegg_counts[z,4]+length(deup);
    kegg_counts[z,5]=kegg_counts[z,5]+length(dedown);
  }
  
  for(kegcat in 1:nrow(keggcat)){
    z=rownames(keggcat)[kegcat];
    deup<-intersect(as.matrix(keggcat[z,3:(2+keggcat[z,1])]),as.list(UP_genes));
    dedown<-intersect(as.matrix(keggcat[z,3:(2+keggcat[z,1])]),as.list(DOWN_genes));
    de<-intersect(as.matrix(keggcat[z,3:(2+keggcat[z,1])]),as.list(DE_genes));
    if((network=="global")||(network=="functional")){
      if(length(de>0)){
        for(f in 1:length(de)){
          Network=rbind(Network,c(z,de[f]));
        }
      }
    }
    keggcat_counts[z,3]=keggcat_counts[z,3]+length(de);
    keggcat_counts[z,4]=keggcat_counts[z,4]+length(deup);
    keggcat_counts[z,5]=keggcat_counts[z,5]+length(dedown);
  }
  
  for(go in 1:nrow(gos)){
    z=rownames(gos)[go];
    deup<-intersect(as.matrix(gos[z,3:(2+gos[z,1])]),as.list(UP_genes));
    dedown<-intersect(as.matrix(gos[z,3:(2+gos[z,1])]),as.list(DOWN_genes));
    de<-intersect(as.matrix(gos[z,3:(2+gos[z,1])]),as.list(DE_genes));
    if((network=="global")||(network=="functional")){
      if(length(de>0)){
        for(f in 1:length(de)){
          Network=rbind(Network,c(z,de[f]));
        }
      }
    }
    go_counts[z,3]=go_counts[z,3]+length(de);
    go_counts[z,4]=go_counts[z,4]+length(deup);
    go_counts[z,5]=go_counts[z,5]+length(dedown);
  }
  
  if(length(which(TF_counts[,3]!=0))==0){
    print("No TF's targets found deregulated")
    ResultsTF=matrix(nrow=1,ncol=1,data="No TF's targets found deregulated");
  } else {
    # 		print(TF_counts[which(TF_counts[,3]!=0),]);
    ResultsTF=matrix(nrow=length(which(TF_counts[,3]!=0)),ncol=3);
    rownames(ResultsTF)=rownames(TF_counts)[which(TF_counts[,3]!=0)];
    colnames(ResultsTF)=c("DE_qvalue","UP_qvalue","DOWN_qvalue");
    ResultsTFtemp=matrix(nrow=length(which(TF_counts[,3]!=0)),ncol=3);
    rownames(ResultsTFtemp)=rownames(TF_counts)[which(TF_counts[,3]!=0)];
    colnames(ResultsTFtemp)=c("DE_pvalue","UP_pvalue","DOWN_pvalue");
    for(i in 1:length(which(TF_counts[,3]!=0))){
      r=which(TF_counts[,3]!=0)[i];
      ResultsTFtemp[i,1]<-phyper(TF_counts[r,4],Number_of_DE_genes,(TF_counts[r,2]-Number_of_DE_genes), TF_counts[r,1],lower.tail=FALSE);
      ResultsTFtemp[i,2]<-phyper(TF_counts[r,5],Number_of_Up_genes,(TF_counts[r,2]-Number_of_Up_genes), TF_counts[r,1],lower.tail=FALSE);
      ResultsTFtemp[i,3]<-phyper(TF_counts[r,6],Number_of_Down_genes,(TF_counts[r,2]-Number_of_Down_genes), TF_counts[r,1],lower.tail=FALSE);
    }
    ResultsTF[,1]=p.adjust(ResultsTFtemp[,1],method="fdr");
    ResultsTF[,2]=p.adjust(ResultsTFtemp[,2],method="fdr");
    ResultsTF[,3]=p.adjust(ResultsTFtemp[,3],method="fdr");
  }
  
  
  
  if(length(which(kegg_counts[,3]!=0))==0){
    print("No Kegg pathway's member found deregulated..Really? Something has gone wrong.");
    Resultskegg=matrix(nrow=1,ncol=1,data="No Kegg pathway's member found deregulated..Really? Something has gone wrong.");
  } else {
    # 		print(kegg_counts[which(kegg_counts[,3]!=0),]);
    Resultskegg=matrix(nrow=length(which(kegg_counts[,3]!=0)),ncol=3);
    rownames(Resultskegg)=rownames(kegg_counts)[which(kegg_counts[,3]!=0)];
    colnames(Resultskegg)=c("DE_qvalue","UP_qvalue","DOWN_qvalue");
    Resultskeggtemp=matrix(nrow=nrow(kegg_counts[which(kegg_counts[,3]!=0),]),ncol=3);
    rownames(Resultskeggtemp)=rownames(kegg_counts)[which(kegg_counts[,3]!=0)];
    colnames(Resultskeggtemp)=c("DE_pvalue","UP_pvalue","DOWN_pvalue");
    for(i in 1:length(which(kegg_counts[,3]!=0))){
      r=which(kegg_counts[,3]!=0)[i];
      Resultskeggtemp[i,1]<-phyper(kegg_counts[r,3],Number_of_DE_genes,(kegg_counts[r,2]-Number_of_DE_genes), kegg_counts[r,1],lower.tail=FALSE);
      Resultskeggtemp[i,2]<-phyper(kegg_counts[r,4],Number_of_Up_genes,(kegg_counts[r,2]-Number_of_Up_genes), kegg_counts[r,1],lower.tail=FALSE);
      Resultskeggtemp[i,3]<-phyper(kegg_counts[r,5],Number_of_Down_genes,(kegg_counts[r,2]-Number_of_Down_genes), kegg_counts[r,1],lower.tail=FALSE);
    }
    Resultskegg[,1]=p.adjust(Resultskeggtemp[,1],method="fdr");
    Resultskegg[,2]=p.adjust(Resultskeggtemp[,2],method="fdr");
    Resultskegg[,3]=p.adjust(Resultskeggtemp[,3],method="fdr");
  }
  
  if(length(which(keggcat_counts[,3]!=0))==0){
    print("No Kegg pathway category's member found deregulated..Really? Something has gone wrong.");
    Resultskeggcat=matrix(nrow=1,ncol=1,data="No Kegg pathway category's member found deregulated..Really? Something has gone wrong.");
  } else {
    # 		print(keggcat_counts[which(keggcat_counts[,3]!=0),]);
    Resultskeggcat=matrix(nrow=length(which(keggcat_counts[,3]!=0)),ncol=3);
    rownames(Resultskeggcat)=rownames(keggcat_counts)[which(keggcat_counts[,3]!=0)];
    colnames(Resultskeggcat)=c("DE_qvalue","UP_qvalue","DOWN_qvalue");
    Resultskeggcattemp=matrix(nrow=nrow(keggcat_counts[which(keggcat_counts[,3]!=0),]),ncol=3);
    rownames(Resultskeggcattemp)=rownames(keggcat_counts)[which(keggcat_counts[,3]!=0)];
    colnames(Resultskeggcattemp)=c("DE_pvalue","UP_pvalue","DOWN_pvalue");
    for(i in 1:length(which(keggcat_counts[,3]!=0))){
      r=which(keggcat_counts[,3]!=0)[i];
      Resultskeggcattemp[i,1]<-phyper(keggcat_counts[r,3],Number_of_DE_genes,(keggcat_counts[r,2]-Number_of_DE_genes), keggcat_counts[r,1],lower.tail=FALSE);
      Resultskeggcattemp[i,2]<-phyper(keggcat_counts[r,4],Number_of_Up_genes,(keggcat_counts[r,2]-Number_of_Up_genes), keggcat_counts[r,1],lower.tail=FALSE);
      Resultskeggcattemp[i,3]<-phyper(keggcat_counts[r,5],Number_of_Down_genes,(keggcat_counts[r,2]-Number_of_Down_genes), keggcat_counts[r,1],lower.tail=FALSE);
    }		
    Resultskeggcat[,1]=p.adjust(Resultskeggcattemp[,1],method="fdr");
    Resultskeggcat[,2]=p.adjust(Resultskeggcattemp[,2],method="fdr");
    Resultskeggcat[,3]=p.adjust(Resultskeggcattemp[,3],method="fdr");
  }
  
  if(length(which(go_counts[,3]!=0))==0){
    print("No GO's member found deregulated..Really? Something has gone wrong.");
    Resultsgo=matrix(nrow=1,ncol=1,data="No GO's member found deregulated..Really? Something has gone wrong.");
  } else {
    # 		print(head(go_counts[which(go_counts[,3]!=0),]));
    Resultsgo=matrix(nrow=length(which(go_counts[,3]!=0)),ncol=3);
    rownames(Resultsgo)=rownames(go_counts)[which(go_counts[,3]!=0)];
    colnames(Resultsgo)=c("DE_qvalue","UP_qvalue","DOWN_qvalue");
    Resultsgotemp=matrix(nrow=nrow(go_counts[which(go_counts[,3]!=0),]),ncol=3);
    rownames(Resultsgotemp)=rownames(go_counts)[which(go_counts[,3]!=0)];
    colnames(Resultsgotemp)=c("DE_pvalue","UP_pvalue","DOWN_pvalue");
    for(i in 1:length(which(go_counts[,3]!=0))){
      r=which(go_counts[,3]!=0)[i];
      Resultsgotemp[i,1]<-phyper(go_counts[r,3],Number_of_DE_genes,(go_counts[r,2]-Number_of_DE_genes), go_counts[r,1],lower.tail=FALSE);
      Resultsgotemp[i,2]<-phyper(go_counts[r,4],Number_of_Up_genes,(go_counts[r,2]-Number_of_Up_genes), go_counts[r,1],lower.tail=FALSE);
      Resultsgotemp[i,3]<-phyper(go_counts[r,5],Number_of_Down_genes,(go_counts[r,2]-Number_of_Down_genes), go_counts[r,1],lower.tail=FALSE);
    }
    Resultsgo[,1]=p.adjust(Resultsgotemp[,1],method="fdr");
    Resultsgo[,2]=p.adjust(Resultsgotemp[,2],method="fdr");
    Resultsgo[,3]=p.adjust(Resultsgotemp[,3],method="fdr");
  }
  
  Network_final=matrix(ncol=2,data="");
  
  if((network=="global")||(network=="functional")){
    for(interaction in 1:nrow(Network)){
      z1=intersect(Network[interaction,1],rownames(Resultskegg));
      z2=intersect(Network[interaction,1],rownames(Resultskeggcat));
      z3=intersect(Network[interaction,1],rownames(Resultsgo));
      if(length(z1>0)){
        if(Resultskegg[z1,1]<=PV_threshold){
          Network_final=rbind(Network_final,Network[interaction,]);
        }
      } else if(length(z2>0)){
        if(Resultskeggcat[z2,1]<=PV_threshold){
          Network_final=rbind(Network_final,Network[interaction,]);
        }
      } else if(length(z3>0)){
        if(Resultsgo[z3,1]<=PV_threshold){
          Network_final=rbind(Network_final,Network[interaction,]);
        }
      } else {
        Network_final=rbind(Network_final,Network[interaction,]);
      }
    }
  } else {
    Network_final=Network;
  }
  
  write.csv(unique(Network_final),file=paste(output,"Network.csv"),quote=F,row.names=F);
  
  
  if(ncol(ResultsTF)>1){
    tempnames=c("Regulator",colnames(ResultsTF));
    ResultsTF=cbind(rownames(ResultsTF),as.data.frame(ResultsTF));
    colnames(ResultsTF)=tempnames;
  }
  
  if(ncol(Resultskegg)>1){
    tempnames=c("KEGG_pathway",colnames(Resultskegg));
    Resultskegg=cbind(rownames(Resultskegg),as.data.frame(Resultskegg));
    colnames(Resultskegg)=tempnames;
  }
  if(ncol(Resultskeggcat)>1){
    tempnames=c("KEGG_pathway_category",colnames(Resultskeggcat));
    Resultskeggcat=cbind(rownames(Resultskeggcat),as.data.frame(Resultskeggcat));
    colnames(Resultskeggcat)=tempnames;
  }
  if(ncol(Resultsgo)>1){
    tempnames=c("GO",colnames(Resultsgo));
    Resultsgo=cbind(rownames(Resultsgo),as.data.frame(Resultsgo));
    colnames(Resultsgo)=tempnames;
  }
  
  
  if(type_of_output=="html"){
    sortable.html.table(as.data.frame(ResultsTF),paste(output,"TF_Enrichment.html"),page.title=paste(output,"TF_Enrichment"));
    
    sortable.html.table(as.data.frame(Resultskegg),paste(output,"KEGG_Enrichment.html"),page.title=paste(output,"KEGG_Enrichment"));
    sortable.html.table(as.data.frame(Resultskeggcat),paste(output,"KEGG_category_Enrichment.html"),page.title=paste(output,"KEGG_category_Enrichment"));
    sortable.html.table(as.data.frame(Resultsgo),paste(output,"GO_Enrichment.html"),page.title=paste(output,"GO_Enrichment"));
  }else if(type_of_output=="csv"){
    write.csv(as.data.frame(ResultsTF),file=paste(output,"TF_Enrichment.csv"),row.names=F);
    
    write.csv(as.data.frame(Resultskegg),file=paste(output,"KEGG_Enrichment.csv"),row.names=F);
    write.csv(as.data.frame(Resultskeggcat),file=paste(output,"KEGG_category_Enrichment.csv"),row.names=F);
    write.csv(as.data.frame(Resultsgo),file=paste(output,"GO_Enrichment.csv"),row.names=F);
  }else {
    print("Please select type of output to be either \"html\" or \"csv\"");
  }
  
  print("Analysis completed!");
  
  
  
}

RNEA(filename = "data_for_markdown/RNEA_E2.txt",species = "Human",FC_threshold = log(2))

RNEA(filename = "data_for_markdown/RNEA_EGF.txt",species = "Human",FC_threshold = log(2))
