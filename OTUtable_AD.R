#####################################################################
#  Jin-Baek Lee [ CURT program ]                                    #
#  Dimensional reduction intestinal microbial data affecting AD  #
#####################################################################

## Set repositories 
setRepositories(ind = 1:8)

## Set working directory
setwd("Your directory path...")
getwd()

## Load library 
library(readxl)
library(tidyverse)
library(datarium)
library(caret) 
library(dplyr)
library(rpart)
library(rpart.plot)
library(sjlabelled)
library(writexl)
library(xlsx)
library(phyloseq)
library(ggfortify)
library(Rtsne)
library(gridExtra)
library(umap)
library(cluster)
library(devtools)
library(ggbiplot)

#install.packages("ggbiplot","vqv")

Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_291')

require("phyloseq")
require(graphics)
packageVersion("phyloseq")  

##View(data)
##View(metadata)

####################################### Load data ########################################
data <- read_excel("sample_bateria.xlsx")
metadata <- read_excel("ERP106411_metadata.xlsx")

dim(data)
glimpse(data)

dim(metadata)
glimpse(metadata)

####################################### Cleansing data ########################################

######### Handling data #########

# Divide microorganism to taxonomy and group by each level
data<-data %>%
  separate(col ="#SampleID",sep=";",into=c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))

data


#### Create each taxonomy table and cleansing ####

# === Genus OTU table === #
Genus_OTUtable <- data%>%group_by(Domain,Kingdom,Phylum,Class,Order,Family,Genus)%>% 
       summarise(ERR2262736=sum(ERR2262736),
                 ERR2263508=sum(ERR2263508),
                 ERR2263509=sum(ERR2263509),
                 ERR2263510=sum(ERR2263510),
                 ERR2263511=sum(ERR2263511),
                 ERR2263512=sum(ERR2263512),
                 ERR2263513=sum(ERR2263513),
                 ERR2263514=sum(ERR2263514),
                 ERR2263515=sum(ERR2263515),
                 ERR2263516=sum(ERR2263516),
                 ERR2263517=sum(ERR2263517),
                 ERR2263518=sum(ERR2263518),
                 ERR2263519=sum(ERR2263519),
                 ERR2263520=sum(ERR2263520),
                 ERR2263521=sum(ERR2263521),
                 ERR2263522=sum(ERR2263522),
                 ERR2263523=sum(ERR2263523),
                 ERR2263524=sum(ERR2263524),
                 ERR2263525=sum(ERR2263525),
                 ERR2263526=sum(ERR2263526),
                 ERR2263527=sum(ERR2263527),
                 ERR2263528=sum(ERR2263528),
                 ERR2263529=sum(ERR2263529),
                 ERR2263530=sum(ERR2263530),
                 ERR2263531=sum(ERR2263531),
                 ERR2263532=sum(ERR2263532),
                 ERR2263533=sum(ERR2263533),
                 ERR2263534=sum(ERR2263534),
                 ERR2263535=sum(ERR2263535),
                 ERR2263536=sum(ERR2263536),
                 ERR2263537=sum(ERR2263537),
                 ERR2263538=sum(ERR2263538),
                 ERR2263539=sum(ERR2263539)
                 )

# Remove taxonomy variables except Genus
Genus_OTUtable <- Genus_OTUtable[,7:length(Genus_OTUtable)]

# Remove unspecified row in Genus variable
Genus_OTUtable <- Genus_OTUtable[!(Genus_OTUtable$Genus=="g__"),]

# Remove NA row in Genus variable
Genus_OTUtable <- Genus_OTUtable%>%na.omit(Genus)

# Clean Genus variable values
Genus_OTUtable <- Genus_OTUtable%>%
  separate(Genus,sep = "__",into = c(NA,"Genus"))

# Transpose Genus OTU table
name <- Genus_OTUtable$Genus
Genus_OTUtable <- as.data.frame(t(Genus_OTUtable[,-1]))
colnames(Genus_OTUtable) <- name

# Create SampleID variable from rownames
Genus_OTUtable<-Genus_OTUtable%>%
  mutate(SampleID=rownames(Genus_OTUtable))

# Reset rownames
rownames(Genus_OTUtable) <- 1:length(rownames(Genus_OTUtable))

dim(Genus_OTUtable)
#View(Genus_OTUtable)

# write excel file
#write.xlsx(Genus_OTUtable,row.names=T,"G:\\내 드라이브\\CURT\\research\\OTUtable\\Genus_OTUtable_own.xlsx")




#================== Using Phyloseq package to compare with taxonomy table which I directly made ====================#

## Comparing Genus_OTUtable(own) & merge_Genus(phyloseq package) ####
otu_mat <- read_excel("sample_bateria.xlsx")
tax_mat <- data

#View(otu_mat)
#View(tax_mat)

## otu_mat
otu_mat<-otu_mat%>%
  mutate(otu=paste0("OTU",1:length(rownames(otu_mat))))

otu_mat <- otu_mat[,-1]

otu_mat <- otu_mat%>%
  tibble::column_to_rownames("otu")

otu_mat <- as.matrix(otu_mat)

otu_mat <- otu_table(otu_mat,taxa_are_rows = TRUE)

## tax_mat
tax_mat <- tax_mat%>%
  separate(col=Domain,sep = "__",into = c(NA,"Domain"))%>%
  separate(col=Kingdom,sep = "__",into = c(NA,"Kingdom"))%>%
  separate(col=Phylum,sep = "__",into = c(NA,"Phylum"))%>%
  separate(col=Class,sep = "__",into = c(NA,"Class"))%>%
  separate(col=Order,sep = "__",into = c(NA,"Order"))%>%
  separate(col=Family,sep = "__",into = c(NA,"Family"))%>%
  separate(col=Genus,sep = "__",into = c(NA,"Genus"))%>%
  separate(col=Species,sep = "__",into = c(NA,"Species"))

tax_mat <- tax_mat%>%
  mutate(otu=paste0("OTU",1:length(rownames(tax_mat))))

tax_mat <- tax_mat[,c(length(tax_mat),1:length(tax_mat)-1)]
tax_mat <- tax_mat[,1:9]

tax_mat <- tax_mat%>%
  tibble::column_to_rownames("otu")

tax_mat <- as.matrix(tax_mat)

tax_mat <- tax_table(tax_mat)


## phyloseq format
carbom <- phyloseq(otu_mat,tax_mat)

taxglom = tax_glom(carbom,"Genus")
rank_names(carbom)

#View(taxglom@otu_table)
#View(taxglom@tax_table)

Genus_id <- c(taxglom@tax_table[,"Genus"])
Genus_taxo <- as.data.frame(taxglom@otu_table@.Data)

merge_Genus <- cbind(Genus_id,Genus_taxo)
#View(merge_Genus)


## Transpose merge_Genus
name <- merge_Genus$Genus_id
merge_Genus <- as.data.frame(t(merge_Genus[,-1]))
colnames(merge_Genus) <- name

#=================================================================================================================#


######### Handling metadata #########

#View(metadata)
metadata <- metadata%>%
  dplyr::select(c(Run,description))

metadata <- metadata%>%
  rename(SampleID=Run)

metadata <- as.data.frame(metadata)

dim(metadata)
metadata$SampleID

## Check if two data.frame can be join
intersect(rownames(Genus_OTUtable),metadata$SampleID)


######### Join OTUtable and metadata #########
Genus_OTUtable <- left_join(metadata,Genus_OTUtable,by="SampleID")


######### Handling each OTUtable ######### 

# AD column : WT = healthy mice = 0 / ADLPART = AD mice = 1
Genus_OTUtable <- Genus_OTUtable%>%
  rename(AD=description)%>%
  mutate(AD=ifelse(Genus_OTUtable$description =='2month_ADLPAPT',1,0))%>%
  mutate(AD=as.factor(AD))

dim(Genus_OTUtable)
###################################### Dimensional reduction ######################################

######### PCA & t-sne #########

# PCA prcomp is based on SVD algorithms / PCA princomp is based on eigen vectors
# m-by-n metrix A 
# m > n : n-by-n   
# m = n : svd(A)
# m < n : m-by-m  m PCs explains 100% variance so, do not need more PCs (think PCs are orthogonal in plot)

## Apply PCA
Genus_OTUtable_pca <- prcomp(Genus_OTUtable[,3:length(Genus_OTUtable)],center=T,scale. = T) 
biplot(Genus_OTUtable_pca,cex=0.8)
## Determine the number of PCs 
summary(Genus_OTUtable_pca) # In Genus, PC25
#View(Genus_OTUtable)
#View(Genus_OTUtable_pca$x)

## Contributions of original variables to the newly created variable PCs = metadata
Genus.metadata <- Genus_OTUtable_pca$rotation
#plot(Genus.metadata)

## PCA result
Genus_OTUtable.dr <-Genus_OTUtable_pca$x[,1:25]

## Add SampleID & AD column from origin OTUtable and set order
Genus_OTUtable.dr<-cbind(Genus_OTUtable.dr,Genus_OTUtable[,1:2])
Genus_OTUtable.dr<-Genus_OTUtable.dr[,c(dim(Genus_OTUtable.dr)[2]-1,1:25,dim(Genus_OTUtable.dr)[2])]

## Final OTUtable

#View(Genus_OTUtable.dr)  ## DR data
#View(Genus_OTUtable)     ## Original data

## Apply t-sne and draw a plot to each dataset

# DR data
tsne.dr <- Rtsne(Genus_OTUtable.dr,perplexity=10)
tsne.dr

tsne.dr.plot <- data.frame(x = tsne.dr$Y[,1], y = tsne.dr$Y[,2], col = Genus_OTUtable.dr$AD)
p1<-ggplot(tsne.dr.plot ) + geom_point(aes(x=x, y=y, color=col)) + labs(title="PCA+tsne")


# Original data 
tsne.origin <- Rtsne(Genus_OTUtable,perplexity=10)
tsne.origin

tsne.origin.plot <- data.frame(x=tsne.origin$Y[,1],y=tsne.origin$Y[,2],col=Genus_OTUtable$AD)
p2<-ggplot(tsne.origin.plot) + geom_point(aes(x=x,y=y,color=col)) + labs(title="Origin+tsne")


#grid.arrange(p1,p2,ncol=2)


######### UMAP #########

# DR data
umap.dr <- umap(Genus_OTUtable.dr[,2:(length(Genus_OTUtable.dr)-1)])
umap.dr.plot <- data.frame(x=umap.dr$layout[,1],y=umap.dr$layout[,2],col=Genus_OTUtable.dr$AD)
up1<-ggplot(umap.dr.plot) +geom_point(aes(x=x,y=y,color=col)) + labs(title="PCA+umap")

# Original data 
umap.origin <- umap(Genus_OTUtable[,3:length(Genus_OTUtable)])
umap.origin.plot <- data.frame(x=umap.origin$layout[,1],y=umap.origin$layout[,2],col=Genus_OTUtable$AD)
up2<-ggplot(umap.origin.plot) + geom_point(aes(x=x,y=y,color=col)) + labs(title="Origin+umap")

grid.arrange(p1,p2,up1,up2,ncol=2,nrow=2)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  Example dataset @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## URL : https://www.kaggle.com/iabhishekofficial/mobile-price-classification?select=train.csv

example.data <- read.csv("train.csv",sep=",")
example.data <- example.data%>%
  mutate(price_range=as.factor(example.data$price_range))

# outcome variable: price_range
#View(example.data)
dim(example.data)
str(example.data)

## PCA data
example_pca <- prcomp(example.data[,1:(length(example.data)-1)],center=T ,scale. =T )  ## 여기에서 PCA와 동시에 정규화 시키면 군집화가 제대로 되지 않음 
summary(example_pca) # select PC16

biplot(example_pca)
example.dr <- example_pca$x[,1:16]

example.dr <- as.data.frame(example.dr)

## proportion of variance in each PCs
index <- seq(from=2,by=3,length.out=16)

lst <- list()
for(i in index){
  lst<-append(lst,summary(example_pca)[["importance"]][i])
}
lst <- unlist(lst)
lst

for (i in 1:dim(example.dr)[1]){
  for( j in 1:dim(example.dr)[2]){
    example.dr[i,j] <- round(example.dr[i,j]*lst[j],4)
  }
}

#View(example.dr)

## pca tsne 
pca.tsne <- Rtsne(example.dr)

pca.tsne.plot <- data.frame(x = pca.tsne$Y[,1], y = pca.tsne$Y[,2], col = example.data$price_range)
p1<-ggplot(pca.tsne.plot) + geom_point(aes(x=x, y=y, color=col)) + labs(title="PCA+tsne")

# Original data 
ori.tsne <- Rtsne(example.data[,1:(length(example.data)-1)],normalize=T)

ori.tsne.plot <- data.frame(x=ori.tsne$Y[,1],y=ori.tsne$Y[,2],col=example.data$price_range)
p2<-ggplot(ori.tsne.plot) + geom_point(aes(x=x,y=y,color=col)) + labs(title="Origin+tsne")


######### UMAP #########

# pca umap
pca.umap <- umap(example.dr,controlscale=T)

pca.umap.plot <- data.frame(x=pca.umap$layout[,1],y=pca.umap$layout[,2],col=example.data$price_range)
up1<-ggplot(pca.umap.plot) +geom_point(aes(x=x,y=y,color=col)) + labs(title="PCA+umap")

# Original data 
ori.umap <- umap(example.data[,1:(length(example.data)-1)],controlscale=T)

ori.umap.plot <- data.frame(x=ori.umap$layout[,1],y=ori.umap$layout[,2],col=example.data$price_range)
up2<-ggplot(ori.umap.plot) + geom_point(aes(x=x,y=y,color=col)) + labs(title="Origin+umap")


grid.arrange(p1,p2,up1,up2,ncol=2,nrow=2)

