#required packages
library(ggplot2)
library(reshape2)
library(fields)

###################
bin_individual_regions <- function(in_data,anno_data,upstream,down,binsize){	
  #find size of all elements in annotation
  anno_data$V1=factor(anno_data$V1)
  size=anno_data$V3 - anno_data$V2
  downstream=(round(max(size),digits=-1)) + down
  #define the number of bins to create based on the largest annotated feature
  num_bin=((upstream + downstream)/binsize)
  
  #initialize output matrix        
  out.matrix=matrix(NA,nrow=nrow(anno_data),ncol=num_bin)
  
  #pull input data based on given upstream and size values. Correct for strandedness
  for(q in 1:length(table(anno_data$V1))){
    in_chr=subset(in_data,as.character(in_data$V1)==names(table(anno_data$V1))[q])
    anno_chr=subset(anno_data,anno_data$V1==names(table(anno_data$V1))[q])
    for(i in 1:nrow(anno_chr)){
      if(anno_chr[i,6]=='+'){
        test=subset(in_chr,in_chr$V2 >= (anno_chr[i,2]-upstream) & in_chr$V2 <= (anno_chr[i,2]+downstream))
        test.bin=stats.bin(test$V2,test$V4,N=num_bin)
        out.line=cbind(matrix(test.bin$centers,ncol=1),test.bin$stats["median",])
        out.matrix[i,1:num_bin]=out.line[1:num_bin,2]
      }
      if(anno_chr[i,6]=='-'){
        test=subset(in_chr,in_chr$V2 >= (anno_chr[i,2]-downstream) & in_chr$V2 <= (anno_chr[i,2]+upstream))
        test.bin=stats.bin(test$V2,test$V4,N=num_bin)
        out.line=cbind(matrix(test.bin$centers,ncol=1),test.bin$stats["median",])
        out.line=out.line[order(-1:-nrow(out.line)),]
        out.matrix[i,1:num_bin]=out.line[1:num_bin,2]
      }
    }
  }
  out.matrix=cbind.data.frame(as.character(anno_data$V4),size,out.matrix)
  colnames(out.matrix)[1:2]=c('ID','size')
  return(out.matrix)
}

#############
omit_distant_bins <- function(binned_data,anno_data,up,down,binsize){
  size=anno_data$V3 - anno_data$V2
  
  filtered.matrix=binned_data
  for(i in 1:nrow(filtered.matrix)){
    if(size[i] != max(size)){
      filtered.matrix[i,(((round(size[i],digits=-1) + up + down) / binsize) + 3):ncol(filtered.matrix)]=NA
    }
  }
  return(filtered.matrix)
}
##############

normalize_by_anno_element <-function(bin_data){
  norm=t(apply(bin_data[,3:ncol(testing.na)],1,function(x){(x / mean(x,na.rm=T))}))
  norm=cbind.data.frame(bin_data[,1:2],norm)
  norm=subset(norm,rowSums(norm[,3:ncol(norm)],na.rm=T)!=0)
  return(norm)
}


test=matrix(seq(1:50),5,10)


##############
create_plot <- function(bin_data,xlabel,ylabel,title){
  melted_data=melt(bin_data,id.vars=c('ID','size'))
  melted_data$value=log(melted_data$value + 1)
  #properly order the ID factors by element size
  melted_data$ID = factor(melted_data$ID,levels=(melted_data$ID)[order(melted_data$size)])
  
  png(paste(title,'.png',sep=''),width=10000,height=10000)
  output_plot=ggplot(melted_data,aes(x=variable,y=ID,fill=value)) + 
    geom_tile() + 
    scale_fill_distiller(type='div',palette = "Spectral",  trans = 'reverse', guide='colourbar') +
    geom_vline(aes(xintercept=100),size=5) +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(title)
  print(output_plot)
  dev.off()
}
##############

#read in annotation (anno) and data of interest (input)
anno=read.delim('TAIR10_primary_genes_simple.bed',head=F)
anno=subset(anno,anno$V1=='Chr1')
input=read.delim('Sample_277_1_WT.plus.real.bed',head=F,nrow=4000000)

a.test=anno[1:40,]

#bin it (the time consuming step)
testing=bin_individual_regions(input,a.test,1000,1000,10)
#filter out the bins beyond the downstream limit for each annotated element
testing.na=omit_distant_bins(testing,a.test,1000,1000,10)
#normalize values by mean of each elements remaining bins
testing.na.normalized=normalize_by_anno_element(testing.na)
#create png heatmap
create_plot(testing.na.normalized,'10bp bins','TAIR10 primary genes','wtplus_norm')


