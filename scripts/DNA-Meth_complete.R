
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: OS X El Capitan 10.11.6
# # Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib


library(readr)
library(multtest)
library(genefilter)

#### write  function for multiple test correction
correctPvalueandReturnAll<-function(tt.pval,method)
{
  mt=mt.rawp2adjp(tt.pval,proc=method)
  adjp=mt$adjp[order(mt$index),]
  return(adjp[,2])
}


##### load the dna.meth-Seq data #####
data.path="gdc_download_20181221_133906.170339"
files <- list.files(path=data.path,recursive=T, pattern = "txt")


# read the first file for the first time
file=files[1]
file.id=strsplit(file,"/")[[1]][1]
#read the first file
temp <- read.table(file.path(data.path,files[1]), sep = "\t", header=T)


#create a storing object dna.meth to save the whole read counts of each file read in an iteration
dna.meth=as.data.frame(temp[,2])
rownames(dna.meth)=temp[,1]
colnames(dna.meth)=c(file.id)


for(i in 2: length(files))
{
  ## refer to the next file (note that we start from index 2, bec we already read the first file)
  file=files[i]
  file.id=strsplit(file,"/")[[1]][1]
  
  # read the next file  
  temp2 <- read.table(file.path(data.path,files[i]), sep = "\t", header=T)

    ## remove the first column, bec we had it already
  temp2=temp2[2]
  colnames(temp2)=c(file.id)
  dna.meth=cbind(dna.meth,temp2)
}

# remove Null rows
dna.meth=dna.meth[complete.cases(dna.meth), ]

#upload pheno file
pheno <- read_delim("gdc_sample_sheet.2018-12-21.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
table(pheno$`Sample Type`)

# rename column names : replace the spaces with dots 
pheno.names=names(pheno)
names(pheno)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))
table(pheno$Sample.Type)

#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids to 
file.ids=colnames(dna.meth)

file.ids.pheno=pheno$File.ID
index.files=match(file.ids,file.ids.pheno)
names(dna.meth)=pheno$Sample.ID[index.files]

#adjust columns if the tumor/normal not beside each others
x = dna.meth[,grep("01A$", colnames(dna.meth))]
y = dna.meth[,grep("11A$", colnames(dna.meth))]
dna.meth2 = cbind(y,x)
table(rownames(dna.meth) == rownames(dna.meth2))

save(dna.meth,pheno,file="DNA-Methy-seq.RDATA")


#### do the DIFF Methylation Analysis
meth=dna.meth2
case.indecies=c(1:3)
ctrl.indecies=c(4:dim(meth)[2])

## calculating LFC
lfc.diff=apply(meth,1, function(x)  mean(x[case.indecies]) -mean(x[ctrl.indecies]))
x=as.data.frame(lfc.diff)

## calculating p values
f=factor( c( rep(1, length(case.indecies)) , rep(2, length(ctrl.indecies)) ))
t.pval=rowttests(as.matrix(meth),f)
t.pval=rowttests(as.matrix(meth),f)$p.value
t.pval.adj=correctPvalueandReturnAll(t.pval,"BH")

res=cbind(lfc.diff,t.pval,t.pval.adj)

#####  selection criteria for identifying dmrs
dmrs.res=res[t.pval<0.05,]  ##### identify dmrs based on the significance level only
dmrs.res=res[t.pval.adj < 0.05,]  ##### identify dmrs based on the significance level only
dmrs.res=res[abs(lfc.diff) > log2(1.5),]  ##### identify dmrs based on the LFC only
dmrs.res=res[abs(lfc.diff) > log2(1.5)   & t.pval <0.05,]  ##### identify dmrs based on both  LFC and the significane level 

dmrs= rownames(dmrs.res)

#expression for dmrs
dmrs.meth3=dna.meth2[rownames(dna.meth2) %in% dmrs, ]
#################################################
#you can make a heatmap and volcano

## don't forget to map to gene symbol before exporting the list of DMRs
temp = temp[temp$Composite.Element.REF %in% dmrs,]
temp = temp[,c("Composite.Element.REF","Gene_Symbol")]
x = temp[!temp$Gene_Symbol == ".",]$Gene_Symbol

 x = temp[!temp$Gene_Symbol == ".",]
 x = x$Gene_Symbol

x = paste0( x, collapse=";")
x = unlist(strsplit(x,";"))
x = x[!duplicated(x)]

#export them for further analysis in DAVID
write.table(x,file = "dmrs_genes.txt",row.names = F,col.names = F,quote = F)

