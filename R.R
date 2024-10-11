library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(data.table)

logFCfilter=1               
adjPvalFilter=0.05          


load("pbmc.rda")


pdf(file="01.featureViolin.pdf", width=14, height=6)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


pdf(file="01.featureCor.pdf",width=12,height=6)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()


top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf",width=10,height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

pdf(file="02.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()


pdf(file="02.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()


pdf(file="02.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()


pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()

ElbowPlot(object = pbmc, ndims = 20)

save(pbmc,file = "pbmc_pca")


pdf(file="03.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    
dev.off()
write.table(pbmc$seurat_clusters,file="03.tsneCluster.txt",quote=F,sep="\t",col.names=F)

 pbmc.markers <- FindAllMarkers(object = pbmc,
                                 only.pos = FALSE,
                                 min.pct = 0.25,
                                 logfc.threshold = logFCfilter)
 
 sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
 write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

 top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

 pdf(file="03.tsneHeatmap.pdf",width=20,height=13)
 DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
 dev.off()



showGenes=c("IL32","CRIP2","ANXA2","VWF")


pdf(file="03.markerViolin.pdf",width=10,height=10)
VlnPlot(object = pbmc, features = showGenes)
dev.off()


pdf(file="03.markerScatter.pdf",width=10,height=10)
FeaturePlot(object = pbmc, features = showGenes, cols = c("green", "red"))
dev.off()


pdf(file="03.markerBubble.pdf",width=12,height=6)
cluster10Marker=showGenes
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()


counts<-pbmc@assays$RNA@counts
clusters<-pbmc@meta.data$seurat_clusters
ann=pbmc@meta.data$orig.ident

 
 ref<-celldex::BlueprintEncodeData()

singler=SingleR(test=counts, ref =ref,
                labels=ref$label.main, clusters = clusters)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
singler2=SingleR(test=counts, ref =ref,
                 labels=ref$label.main)
cellAnn=as.data.frame(singler2)
cellAnn=cbind(id=row.names(cellAnn), cellAnn)
cellAnn=cellAnn[,c("id", "labels")]
write.table(cellAnn, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)



newLabels=singler$labels
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, newLabels)
pdf(file="04.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    
dev.off()


save(pbmc,file = "pbmc1.rda")

load("pbmc1.rda")


pbmc.markers=FindAllMarkers(object = pbmc,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.cellMarkers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.cellMarkers,file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)



library(SeuratData)
library(patchwork)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(Seurat)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)



getScatterplot <- function(object, gene1, gene2, cor.method = "pearson",
                           jitter.num = 0.15, pos = TRUE){
  if (!gene1 %in% rownames(object)) {
    print("gene1 was not found")
    if (!gene2 %in% rownames(object)) {
      print("gene2 was not found")
    }
  }else{
    exp.mat <- GetAssayData(object = object, assay = "RNA") %>% .[c(gene1,gene2),] %>% 
      as.matrix() %>% t() %>% as.data.frame()
    if (pos) {
      if(nrow(exp.mat[which(exp.mat[,1] > 0 & exp.mat[,2] > 0),]) > (nrow(exp.mat)*0.01)){
        exp.mat <- exp.mat[which(exp.mat[,1] > 0 & exp.mat[,2] > 0),]
      }else{
        exp.mat <- exp.mat[which(exp.mat[,1] > 0 | exp.mat[,2] > 0),]
      }
    }
    colnames(exp.mat) <- c("Var1", "Var2")
    plots <- ggplot(data=exp.mat, mapping = aes_string(x = "Var1", y = "Var2")) + 
      geom_smooth(method = 'lm', se = T, color='red', size=1) +
      stat_cor(method = cor.method)+ labs(x=gene1, y=gene2) +
      geom_jitter(width = jitter.num, height = jitter.num, color = "black", size = 1, alpha=1)+
      theme_bw()+
      theme(panel.grid=element_blank(),
            legend.text=element_text(colour= 'black',size=10),
            axis.text= element_text(colour= 'black',size=10),
            axis.line= element_line(colour= 'black'),
            panel.border = element_rect(size = 1, linetype = "solid", colour = "black"),
            panel.background=element_rect(fill="white"))
    return(plots)
  }
}


showGenes <- c("GLS","DLD","LIPT1")
geneCard <- intersect(c("CDKN2A","FDX1","DLAT","LIAS","MTF1","PDHA1","PDHB"),rownames(pbmc))

for (j in geneCard) {
  for (i in showGenes) {
    
    p1 <- FeaturePlot(pbmc, features = c(i, j), 
                      blend = TRUE, cols = c("gray80","red", "green"), 
                      pt.size = 0.5, raster = F) +  
      theme(aspect.ratio = 1)
    p2 <- getScatterplot(pbmc, gene1 = j, gene2 = i, 
                         jitter.num = 0.15, pos = TRUE) +
      theme(aspect.ratio = 1)
    p <- CombinePlots(plots = list(p1, p2), ncol = 2, rel_widths = c(4, 1))
    pdf(file=paste0(j," ~ ",i,".pdf"),width=15,height=4)
    print(p)
    dev.off()
  }
}


library(glmnet)
library(survival)
library(pheatmap)
library(gplots)
library(survcomp)
library(survivalROC)
library(pROC)
library(ggplot2)

display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
  
} 


jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")









workdir <- "F:\\工作文档\\非肿瘤\\胆汁淤积\\4.Lasso建模\\R"; setwd(workdir)




tpms <- read.table("ARGexp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpms[is.na(tpms)] <- 0
colnames(tpms) <- gsub("-","_",colnames(tpms))
tum.sam <- rownames(tpms)

geo1.tpms <- read.csv("1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
geo1.tpms[is.na(geo1.tpms)] <- 0
colnames(geo1.tpms) <- gsub("-","_",colnames(geo1.tpms))


comgene <- intersect(colnames(tpms),colnames(geo1.tpms))

length(comgene)
tpms <- tpms[,comgene]
geo1.tpms <- geo1.tpms[,comgene]







outTab <- NULL
for (seed in 1:100) {
  
  
  risk <- NULL
  
  set.seed(seed = seed)
  tmp <- tpms[tum.sam,]
  colnames(tmp) <- make.names(colnames(tmp))
  
  
  cvfit = cv.glmnet(as.matrix(tmp[,-1]),
                    tmp$Tissue,
                    family = "gaussian",
                    alpha = 1,
                    nfold = 5) 
  myCoefs <- coef(cvfit, s='lambda.min');
  fea <- rownames(coef(cvfit, s = 'lambda.min'))[coef(cvfit, s = 'lambda.min')[,1]!= 0]
  
  lasso_fea <- fea
  lasso_coef <- myCoefs@x; names(lasso_coef) <- lasso_fea
  
  lasso_coef.hr <- data.frame(gene = names(lasso_coef),
                              coef = lasso_coef,
                              
                              
                              
                              stringsAsFactors = F)
  
  lasso_coef.hr <- lasso_coef.hr[intersect(comgene,names(lasso_coef)),]
  lasso_coef.hr <- lasso_coef.hr[order(lasso_coef.hr$coef,decreasing = F),]
  
  
  if(nrow(lasso_coef.hr) > 1) {
    
    tmp <- as.matrix(tpms[tum.sam,rownames(lasso_coef.hr)])
    risk.score <- apply(tmp,1,function(x) {x %*% lasso_coef.hr$coef})
    risk.score.train <- risk.score
    
    tmp <- tpms[names(risk.score),1:2]
    
    tmp$risk.score <- as.numeric(risk.score)
    tmp$RiskGroup <- ifelse(tmp$risk.score > median(risk.score) ,"HRisk","LRisk")
    risk <- rbind.data.frame(risk,
                             data.frame(samID = tum.sam,
                                        riskscore = tmp$risk.score,
                                        riskgroup = tmp$RiskGroup,
                                        cohort = "GEO training",
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
    
    fit <- glm(Tissue ~ risk.score, data=tmp, family = "gaussian")
    predicted <- predict(fit, tmp, type = "response")
    rocobj <- roc(tmp$Tissue,tmp$risk.score)
    auc <- round(auc(tmp$Tissue,tmp$risk.score), 4)
    
    
    
    
    tmp.validate <- as.matrix(geo1.tpms[,rownames(lasso_coef.hr)])
    risk.score <- apply(tmp.validate,1,function(x) {x %*% lasso_coef.hr$coef})
    risk.score.validate <- risk.score
    
    tmp.validate <- geo1.tpms[names(risk.score),1:2]
    tmp.validate$risk.score <- as.numeric(risk.score)
    tmp.validate$RiskGroup <- ifelse(tmp.validate$risk.score > median(risk.score) ,"HRisk","LRisk")
    risk <- rbind.data.frame(risk,
                             data.frame(samID = rownames(geo1.tpms),
                                        riskscore = tmp.validate$risk.score,
                                        riskgroup = tmp.validate$RiskGroup,
                                        cohort = "GEO test",
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
    
    predicted.test <- predict(fit, tmp.validate, type = "response")
    rocobj.test <- roc(tmp.validate$Tissue,tmp.validate$risk.score)
    auc.test <- round(auc(tmp.validate$Tissue,tmp.validate$risk.score), 4)
    
    
    if(auc > 0.5 & auc.test > 0.5) {
      cat(paste0("seed=",seed,"; auc.train=",auc,"; auc.test=",auc.test,"\n"))
      cat("\n")
      outTab <- rbind.data.frame(outTab,data.frame(seed=seed,
                                                   
                                                   modelgene.num=nrow(lasso_coef.hr),
                                                   
                                                   auc.rain=auc,
                                                   auc.test=auc.test,
                                                   modelgene=paste(gsub("-","_",rownames(lasso_coef.hr)),collapse = ","),
                                                   
                                                   stringsAsFactors = F),
                                 stringsAsFactors = F)
      
      p1 <- ggroc(rocobj,color = "red",linetype = 1, size = 1, alpha =1,legacy.axes = T) +
        geom_abline(intercept=0,slope=1,color="grey",size=1,linetype=1) +
        labs(x="Specificity",
             y="Sensivity",
             title = "GEO train") +
        annotate("text",x=.8,y=.1,label=paste0("AUC = ",auc),
                 size =5,family="serif") +
        coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
        theme_bw() +
        theme(panel.background = element_rect(fill='transparent'),
              axis.ticks.length = unit(0.4,"lines"),
              axis.ticks = element_line(color='black'),
              axis.line = element_line(size=.5, colour='black'),
              axis.title = element_text(colour='black',size=12,face="bold"),
              axis.text = element_text(colour='black',size=10,face="bold"),
              text = element_text(colour='black',size=8,family="serif"))
      p2 <- ggroc(rocobj.test,color = "red",linetype = 1, size = 1, alpha =1,legacy.axes = T) +
        geom_abline(intercept=0,slope=1,color="grey",size=1,linetype=1) +
        labs(x="Specificity",
             y="Sensivity",
             title = "GEO test") +
        annotate("text",x=.8,y=.1,label=paste0("AUC = ",auc.test),
                 size =5,family="serif") +
        coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
        theme_bw() +
        theme(panel.background = element_rect(fill='transparent'),
              axis.ticks.length = unit(0.4,"lines"),
              axis.ticks = element_line(color='black'),
              axis.line = element_line(size=.5, colour='black'),
              axis.title = element_text(colour='black',size=12,face="bold"),
              axis.text = element_text(colour='black',size=10,face="bold"),
              text = element_text(colour='black',size=8,family="serif"))
      pdf(file = paste0("survivalROC for training dataset",seed,".pdf"),width = 4,height = 4)
      print(p1)
      dev.off()
      pdf(file = paste0("survivalROC for validation dataset",seed,".pdf"),width = 4,height = 4)
      print(p2)
      dev.off()
    }
  }
  write.table(outTab,paste0("outTab",seed,".txt"),sep = "\t",row.names = F,quote = F)
}



write.table(outTab,"outTab.txt",sep = "\t",row.names = F,quote = F) 
write.table(risk,"risk.txt",sep = "\t",row.names = F,quote = F)      






darkred   <- "#F2042C"
darkblue   <- "#21498D"

cutoff <- 0.1
lasso_coef.hr$gene <- gsub("_","-",lasso_coef.hr$gene)
lasso_coef.hr$group <- as.character(cut(lasso_coef.hr$coef, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c("#21498D","#EABF00","#21498D")))
pdf("lasso_coef_hr.pdf",width = 3.5,height = 2.5)
par(bty="n", mgp = c(1.7,.33,0),mar=c(2.5,2.7,1,1)+.1, las=1, tcl=-.25,xpd = T)
a <- barplot(lasso_coef.hr$coef,col = lasso_coef.hr$group,border = NA,
             horiz = T,xlim = c(-1,1),add=F,xaxt = "n")
axis(side = 1, at = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),
     labels = c("-1","-.8","-.6","-.4","-.2","0",".2",".4",".6",".8","1"))



for (i in 1:nrow(lasso_coef.hr)) {
  text(y = a[,1][i],x = ifelse(lasso_coef.hr$coef[i] > 0,-0.0001,0.0001),pos = ifelse(lasso_coef.hr$coef[i] > 0,2,4),labels = lasso_coef.hr$gene[i],adj = ifelse(lasso_coef.hr$coef[i]>0,0,1))
}

points(0.6,1,pch = 19, cex = 1.5)

text(0.6,1,"Coefficient",pos = 4)
invisible(dev.off())
write.table(lasso_coef.hr[,1:2], "lasso coefficient.txt",sep = "\t",row.names = F,col.names = T,quote = F)


pdf("lasso.pdf",width = 4.5,height = 4)
par(bty="o", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25,xpd = F)
plot(cvfit$glmnet.fit, "lambda", label=F)
abline(v=log(cvfit$lambda.min),lwd=2,col="grey60",lty=4)
invisible(dev.off())

pdf("cvfit.pdf",width = 4.5,height = 4)
par(bty="o", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25)
plot(cvfit)
abline(h=min(cvfit$cvm),lwd=2,col="black",lty=4)
points(log(cvfit$lambda.min),min(cvfit$cvm),pch=18,cex=2,col="black")
points(log(cvfit$lambda.min),min(cvfit$cvm),pch=18,cex=1.5,col="#008B8A")
invisible(dev.off())

save.image("Heng.RData")



input="CIBERSORT-Results.txt"
outpdf="barplot.pdf"

data <- read.table(input,header=T,sep="\t",check.names=F,row.names=1)
data=t(data)
col=rainbow(nrow(data),s=0.7,v=0.7)

pdf(outpdf,height=10,width=15)
par(las=1,mar=c(8,4,4,15))
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=65,xpd=T);text(a1,-0.02,colnames(data),adj=1.1,cex=1.1);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])

legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()


rt=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)

library(corrplot)
pdf("corHeatmap.pdf",height=13,width=13)              
corrplot(corr=cor(rt),
         method = "color",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white", "red"))(50),
)
dev.off()


normal=7                                                          
tumor=3                                                         

rt=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)   

pdf("vioplot.pdf",height=8,width=15)              
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,64),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)

for(i in 1:ncol(rt)){
  normalData=rt[1:normal,i]
  tumorData=rt[(normal+1):(normal+tumor),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p < 0.05, ifelse(p < 0.01,"**","*"), ""),cex = 0.8)
  
  
}
dev.off()




library(ggplot2)
library(ggpubr)
library(SimDesign)
library(cowplot)
library(dplyr)
library(GSVA)
library(limma)
library(stringr)


jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")


expr <- read.table("symbol.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- read.table("genelist.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = F)




ciber <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)






for (i in gene$V1) {
  message(paste0("analysis of ",i," starts..."))
  subexpr <- as.numeric(expr[i,])
  names(subexpr) <- colnames(expr)
  
  lsam <- names(subexpr[subexpr < median(subexpr)])
  hsam <- names(subexpr[subexpr >= median(subexpr)])

  
  dat <- as.numeric(expr[i,]); names(dat) <- colnames(expr)
  comsam <- intersect(names(dat), rownames(ciber))
  tmp1 <- dat[comsam]
  tmp2 <- ciber[comsam,]
  
  var <- colnames(ciber)
  data <- data.frame(var)
  for (j in 1:length(var)){
    test <- cor.test(as.numeric(tmp2[,j]),tmp1,method = "spearman") 
    data[j,2] <- test$estimate                                            
    data[j,3] <- test$p.value
  }
  names(data) <- c("symbol","correlation","pvalue")
  data <- as.data.frame(na.omit(data))
  data %>% 
    filter(pvalue <0.05) %>%  
    ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
    geom_segment(aes(xend=0,yend=symbol)) +
    geom_point(aes(col=pvalue,size=abs(correlation))) +
    scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
    scale_size_continuous(range =c(2,8))  +
    theme_minimal() +
    ylab(NULL)
  ggsave(paste0("correlation between cibersort and expression of ", i,".pdf"),width = 8,height = 6)
}


library(ggplot2)
library(stringr)


a_1 <- read.table("symbol.txt",header = T,row.names = 1,sep = "\t", quote = "",fill = T,check.names=F)
dim(a_1)
head(a_1[,1:3])
a_2 <- as.data.frame(t(a_1))
dim(a_2)
head(a_2[,1:3])

a_3 <- a_1
a_3$Id <- rownames(a_3)
dim(a_3)
head(a_3[,1:3])

b_1 <- read.table("Immunomodulator_and_chemokines.txt",header = T,sep = "\t", quote = "",fill = T)
dim(b_1)
head(b_1)

b_2 <- b_1[b_1$type == "chemokine",]
dim(b_2)
head(b_2)

data1 <- dplyr::inner_join(b_2,a_3,by="Id")
dim(data1)
head(data1[,1:6])
data2 <- a_2[,c("IL32","CRIP2","ANXA2","VWF",data1$Id)]
dim(data2)
head(data2[,1:5])



library(Hmisc)

CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor) 
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] )
}

res <- rcorr(as.matrix(data2),type = "pearson")
result_1 <- CorMatrix(res$r, res$P)
head(result_1)

dim(result_1)     
result_2 <- result_1[result_1$row == "IL32" |result_1$row == "CRIP2"|result_1$row == "ANXA2"|result_1$row == "VWF",]
dim(result_2)
head(b_2)
b_2$column <- b_2$Id
head(b_2)
result_3 <- dplyr::inner_join(result_2,b_2,by="column")
dim(result_3)
result1 <- result_3[,1:4]
head(result1)
dim(result1)
result1$Regulation <- result1$cor
result1[,5][result1[,5] > 0] <- c("postive")
result1[,5][result1[,5] < 0] <- c("negative")
head(result1)
colnames(result1) <- c("gene", "immuneGene", "cor", "pvalue", "Regulation")

write.table(result1,file="chemokine.xls",sep="\t",quote=F,col.names=T,row.names = F)



a1 <- read.table("chemokine.xls",header = T,sep = "\t", quote = "",fill = T)
head(a1)

data2 <- a1
library(ggpubr)
data2$pvalue <- ifelse(data2$pvalue < 0.05,
                       ifelse(data2$pvalue < 0.01,"**","*"),
                       "")
data2$pvalue[1:20]
data2$type <- data2$cor


summary(data2)
data3 <- data2[order(data2$immuneGene,data2$cor),]
head(data3)
dim(data3)
data4 <- data3[data3$pvalue < 0.05,]
dim(data4)
summary(data4)

p <- ggplot(data4,aes(x=gene,y=immuneGene)) +
  geom_point(aes(colour = cor, size=pvalue)) +
  labs(x="",y="chemokine")
p <- p + scale_colour_gradient2(low = "blue", high = "red", mid = "white",
                                midpoint = 0, limit = c(-1, 1), space = "Lab",
                                name="Pearson\nCorrelation")
p <- p + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  theme(axis.text=element_text(size = 15)) +
  theme(axis.text.x=element_text(colour = "black",angle=0,hjust=0.5,size = 15)) +
  theme(axis.text.y=element_text(colour = "black", vjust=0,size = 15)) +
  theme(axis.title =element_text(size = 20)) +
  theme(text = element_text(size = 15))

p+rotate_x_text(45)
ggsave("chemokine.pdf")



library(cowplot)
library(clusterProfiler)
library(enrichplot)
library(plyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 





selectedGeneID <- c("IL32")                             


mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")



gsym.fc <- read.table("input.txt", header = T)
dim(gsym.fc)





gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")





gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)



gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]



id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID






kk <- gseKEGG(id.fc, organism = "hsa")





kk.gsym <- setReadable(kk, 'org.Hs.eg.db', 
                       'ENTREZID')


sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T),]




write.csv(sortkk,"gsea_output.csv", quote = F, row.names = F)               



geneSetID <- c("hsa05330", "hsa04145", "hsa04064")





for (i in geneSetID) {
  gseaplot(kk, i)
  myGeneList <- enrichplot:::gsInfo(kk, i)
  row.names(myGeneList) <- gsym.fc$gsym
  myGeneList$id <- gsym.fc$ENTREZID 
  write.csv(myGeneList, paste0("gsea_genelist_", i, "_group1.csv"))
}

x <- kk
geneList <- position <- NULL 


gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL,3)                                 


p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
  scale_color_manual(values = mycol) +
  
  
  geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) + 
  ylab("Enrichment\n Score") +
  
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  
  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +
  
  theme(axis.text.y=element_text(size = 12, face = "bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))




rel_heights <- c(1.5, .5, 1.5) 

i <- 0
for (term in unique(gsdata$Description)) {
  idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
  gsdata[idx, "ymin"] <- i
  gsdata[idx, "ymax"] <- i + 1
  i <- i + 1
}


p2 <- ggplot(gsdata, aes_(x = ~x)) +
  geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
  xlab(NULL) + ylab(NULL) + 
  scale_color_manual(values = mycol) + 
  
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  
  theme(legend.position = "none",
        plot.margin = margin(t=-.1, b=0,unit="cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank()) +
  
  scale_y_continuous(expand=c(0,0))




df2 <- p.res$data
df2$y <- p.res$data$geneList[df2$x]
df2$gsym <- p.res$data$gsym[df2$x]



selectgenes <- data.frame(gsym = selectedGeneID)
selectgenes <- merge(selectgenes, df2, by = "gsym")
selectgenes <- selectgenes[selectgenes$position == 1,]
head(selectgenes)


p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) + 
  geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), 
               color = "grey") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycol, guide=FALSE) + 
  scale_color_manual(values = mycol, guide=FALSE) + 
  
  
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + 
  ylab("Ranked list\n metric") +
  xlab("Rank in ordered dataset") +
  
  theme_bw() +
  theme(axis.text.y=element_text(size = 12, face = "bold"),
        panel.grid = element_blank()) +
  
  
geom_text_repel(data = selectgenes, 
                show.legend = FALSE, 
                direction = "x", 
                ylim = c(2, NA), 
                angle = 90, 
                size = 2.5, box.padding = unit(0.35, "lines"), 
                point.padding = unit(0.3, "lines")) +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))




plotlist <- list(p.res, p2, p.pos)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(size = 12, face = "bold"))

plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)



ggsave("GSEA_multi_pathways.pdf", width=6, height=5)





library(clusterProfiler)
library(GOplot)
library(tidyverse)
library(data.table)
library(ggraph)
library(tidygraph)
library(dplyr)






sortkk <- kk.gsym[kk.gsym@result$Description %like% "Allograft rejection" | 
                    kk.gsym@result$Description %like% "Phagosome" | 
                    kk.gsym@result$Description %like% "NF-kappa B signaling pathway",]  


go <- data.frame(Category = "KEGG",
                 ID = sortkk$ID,
                 Term = sortkk$Description, 
                 Genes = gsub("/", ", ", sortkk$core_enrichment), 
                 adj_pval = sortkk$p.adjust)


genelist <- data.frame(ID = gsym.fc.id$SYMBOL, logFC = gsym.fc.id$logFC)


circ <- circle_dat(go, genelist)
head(circ)





write.csv(circ[,c(3,5,6)],"very_easy_input.csv", quote = F, row.names = F)



df <- read.csv("very_easy_input.csv")
head(df)
source(file = "gather_graph_node.R")
source(file = "gather_graph_edge.R")
nodes <- gather_graph_node(df, index = c("term", "genes"), value = "logFC", root="all")
edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
nodes <- nodes %>% mutate_at(c("node.level","node.branch"),as.character)
head(nodes, 10)




graph <- tbl_graph(nodes, edges)


gc <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
  
  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 1/3,edge_width=1) + 
  geom_node_point(aes(size = node.size, 
                      color = node.branch,
                      filter=node.level!="all"), alpha = 1/3) + 
  scale_size(range = c(0.5,20)) + 
  theme(legend.position = "none") + 
  
  
  
  scale_edge_color_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") +
  
  
  geom_node_text(
    aes(
      x = 1.048 * x, 
      y = 1.048 * y, 
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf,
      color = node.branch
    ),
    size = 3, hjust = 'outward') +
  
  
  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all"),
        color = node.branch),
    fontface="bold",
    size=6,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) +
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) 

gc


ggsave("ccgraph_color.pdf", width = 14, height = 14)


gc1 <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
  
  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 0.5, 
                     edge_width=2.5) + 
  scale_edge_color_manual(values = c("#61C3ED","red","purple","darkgreen")) + #?Զ?????ɫ
  
  
  geom_node_point(aes(size = node.size,
                      filter=node.level!="all"), 
                  
                  color = "#61C3ED") + #ͳһΪ??��ɫ
  scale_size(range = c(0.5,20)) + 
  theme(legend.position = "none") + 
  
  
  geom_node_text(
    aes(
      x = 1.05 * x, 
      y = 1.05 * y, 
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf
    ),
    color="black", 
    size = 6, hjust = 'outward') +
  
  
  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all")
    ),
    color="black", 
    fontface="bold",
    size=6,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) + 
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) 

gc1


ggsave("ccgraph.pdf",width = 14,height = 14)



library(ggplot2)
library(ggpubr)
library(SimDesign)
library(cowplot)
library(dplyr)
library(GSVA)
library(limma)
library(stringr)


jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")


expr <- read.table("symbol.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- read.table("genelist.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = F)




ciber <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
load("hallmark.gs.RData")










colnames(expr) <- substr(colnames(expr), 1, 10)


gsva_es <- gsva(as.matrix(log2(expr + 1)), gs)
cutoff <- 1

for (i in gene$V1) {
  message(paste0("analysis of ",i," starts..."))
  subexpr <- as.numeric(expr[i,])
  names(subexpr) <- colnames(expr)
  
  tmp<-subexpr
  median.cut <- median(tmp)
  hsam <- names(tmp[tmp > median.cut])
  lsam <- names(tmp[tmp <= median.cut])
  
  
  dat <- as.numeric((expr[i,])); names(dat) <- colnames(expr)
  comsam <- intersect(names(dat), rownames(ciber))
  tmp1 <- dat[comsam]
  tmp2 <- ciber[comsam,]
  
  var <- setdiff(colnames(ciber),"CancerType")
  data <- data.frame(var)
  for (j in 1:length(var)){
    test <- cor.test(as.numeric(tmp2[,j]),tmp1,method = "pearson") 
    data[j,2] <- test$estimate
    data[j,3] <- test$p.value
  }
  names(data) <- c("symbol","correlation","pvalue")
  data <- as.data.frame(na.omit(data))
  data %>%
    filter(pvalue <0.05) %>%  
    ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
    geom_segment(aes(xend=0,yend=symbol)) +
    geom_point(aes(col=pvalue,size=abs(correlation))) +
    scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
    scale_size_continuous(range =c(2,8))  +
    theme_minimal() +
    ylab(NULL)
  ggsave(paste0("correlation between cibersort and expression of ", i,".pdf"),width = 8,height = 6)
  
  
  group_list <- data.frame(sample = c(lsam,hsam), group = c(rep("LExp", length(lsam)), rep("HExp", length(hsam))))
  design <- model.matrix(~ 0 + factor(group_list$group))
  colnames(design) <- levels(factor(group_list$group))
  rownames(design) <- colnames(gsva_es)
  contrast.matrix <- makeContrasts(HExp-LExp, levels = design)
  
  fit <- lmFit(gsva_es[,group_list$sample], design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  
  pathway <- str_replace(row.names(x), "HALLMARK_", "")
  df <- data.frame(ID = pathway, score = x$t)
  df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))
  df$ID <- gsub("_"," ",paste0(substr(df$ID,1,1), tolower(substr(df$ID,2,nchar(df$ID)))))
  sortdf <- df[order(df$score),]
  sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
  
  ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') +
    coord_flip() +
    scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) +
    geom_hline(yintercept = c(-cutoff,cutoff),
               color="white",
               linetype = 2,
               size = 0.3) +
    geom_text(data = subset(df, score < 0),
              aes(x=ID, y= 0, label= paste0(" ", ID), color = group),
              size = 3,
              hjust = "inward" ) +
    geom_text(data = subset(df, score > 0),
              aes(x=ID, y= -0.1, label=ID, color = group),
              size = 3, hjust = "outward") +
    scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
    
    xlab("") +ylab(paste0("t value of GSVA score\n HExp vs LExp group of ",i)) +
    theme_bw() +
    theme(panel.grid =element_blank()) +
    theme(panel.border = element_rect(size = 0.6)) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  ggsave(paste0("GSVA plot of ", i,".pdf"),width = 8,height = 7)
}


library(AUCell)
library(RcisTarget)
library(doMC)
library(doRNG)
library(DT)
library(visNetwork)





data("motifAnnotations_hgnc_v9")





gene <- read.table("keygene.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
geneList1 <- rownames(gene)
head(geneList1)



geneLists <- list(key_gene=geneList1)


data(motifAnnotations_hgnc)





library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
data(hg19_500bpUpstream_motifRanking_cispbOnly)
motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly
motifRankings


motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_hgnc_v9)
head(motifEnrichmentTable_wGenes)


motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)

resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]


datatable(resultsSubset[,-c("rankAtMax","TF_lowConf"), with=FALSE],       
          escape = FALSE, 
          filter="top", options=list(pageLength=5))


write.table(data.frame('ID'=row.names(motifEnrichmentTable_wGenes),motifEnrichmentTable_wGenes),file='onco_matrix1.txt',sep='\t',quote=F,row.names = F)


motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)
(motifs_AUC )

par(mfrow = c(1,2))

pdf(file="motif enrichment.pdf",width=8,height=8)
for(i in names(geneLists)){
  auc <- getAUC(motifs_AUC)[i,]
  hist(auc, main=i, xlab="AUC histogram",
       breaks=100, col="#ff000050", border="darkred")
  nes3 <- (3*sd(auc)) + mean(auc)
  abline(v=nes3, col="red")
}
dev.off()


motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                           motifAnnot=motifAnnotations_hgnc_v9,
                                           highlightTFs=list(lasso_gene="CEBPB"))                    
head(motifEnrichmentTable[,-"TF_lowConf", with=FALSE])

motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                   rankings=motifRankings, 
                                                   geneSets=geneLists)

motifEnrichmentTable_wGenes[1:4,]

geneSetName <- names(geneLists)[1]

selectedMotifs <- c("cisbp__M1854","cisbp__M4879","cisbp__M0624")     



pdf(file="motif enrichment best.pdf",width=8,height=8)
par(mfrow=c(2,2))
getSignificantGenes(geneLists[[geneSetName]], 
                    motifRankings,
                    signifRankingNames=selectedMotifs,
                    plotCurve=TRUE, maxRank=5000, genesFormat="none",
                    method="aprox")
dev.off()



library(ggpubr)

pFilter=0.99

rt=read.table("ARGexp.txt",sep="\t",header=T,row.names=1,check.names=F)    




data=rt




Type=read.table("cluster.Immunity.txt",sep="\t",check.names=F,row.names=1,header=F)
Type=Type[row.names(data),]
colnames(Type)=c("cluster","Subtype")

outTab=data.frame()
data=cbind(data,Type)
for(i in colnames(data[,1:(ncol(data)-2)])){
  rt1=data[,c(i,"Subtype")]
  colnames(rt1)=c("expression","Subtype")
  ksTest<-kruskal.test(expression ~ Subtype, data = rt1)
  pValue=ksTest$p.value
  if(pValue<pFilter){
    outTab=rbind(outTab,cbind(rt1,gene=i))
    print(pValue)
  }
}
write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)

data=read.table("data.txt",sep="\t",header=T,check.names=F)       
data$Subtype=factor(data$Subtype, levels=c("Control","Disease"))
p=ggboxplot(data, x="gene", y="expression",color = "grey",fill = "Subtype",
            ylab="Expression",
            xlab="",
            palette =c("skyblue","pink") ) 
p=p+rotate_x_text(45)
p
pdf(file="boxplot.pdf",width=12,height=4)                          

p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",method="anova")
dev.off()

library(dplyr)
library(ggplot2)
data %>% 
  filter(Subtype %in% c("Control","Disease")) %>% 
  ggplot(aes(x= gene, y= expression, fill = Subtype, color = Subtype))+
  geom_boxplot(alpha=0.3)+
  scale_fill_manual(name= "Subtype", values = c("deepskyblue", "hotpink"))+
  scale_color_manual(name = "Subtype", values = c("dodgerblue", "plum3"))+
  theme_bw()+labs(x="", y="Expression")+
  theme(axis.text.x = element_text( vjust = 1,size = 12, hjust = 1,colour = "black"),legend.position="top")+
  rotate_x_text(45)+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",method="anova")


library(ggplot2)
library(stringr)




















a_1 <- read.table("symbol.txt",header = T,row.names = 1,sep = "\t", quote = "",fill = T,check.names=F)
a_1=log2(a_1+1)
dim(a_1)
head(a_1[,1:3])
a_2 <- as.data.frame(t(a_1))
dim(a_2)
head(a_2[,1:3])

a_3 <- a_1
a_3$Id <- rownames(a_3)
dim(a_3)
head(a_3[,1:3])

b_1 <- read.table("111.txt",header = T,sep = "\t", quote = "",fill = T)
dim(b_1)
head(b_1)

b_2 <- b_1[b_1$type == "Disease",]
dim(b_2)
head(b_2)

data1 <- dplyr::inner_join(b_2,a_3,by="Id")
dim(data1)
head(data1[,1:6])
data2 <- a_2[,c("IL32", "CRIP2","ANXA2","VWF",data1$Id)]
dim(data2)
head(data2[,1:5])



library(Hmisc)

CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor) 
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] )
}

res <- rcorr(as.matrix(data2),type = "pearson")
result_1 <- CorMatrix(res$r, res$P)
head(result_1)

dim(result_1)     
result_2 <- result_1[result_1$row == "IL32"  |result_1$row == "CRIP2" |result_1$row == "ANXA2"|result_1$row == "VWF",]
dim(result_2)
head(b_2)
b_2$column <- b_2$Id
head(b_2)
result_3 <- dplyr::inner_join(result_2,b_2,by="column")
dim(result_3)
result1 <- result_3[,1:4]
head(result1)
dim(result1)
result1$Regulation <- result1$cor
result1[,5][result1[,5] > 0] <- c("postive")
result1[,5][result1[,5] < 0] <- c("negative")
head(result1)
colnames(result1) <- c("gene", "immuneGene", "cor", "pvalue", "Regulation")

write.table(result1,file="Cholestasis.xls",sep="\t",quote=F,col.names=T,row.names = F)



a1 <- read.table("Cholestasis.xls",header = T,sep = "\t", quote = "",fill = T)
head(a1)

data2 <- a1
library(ggpubr)
data2$pvalue <- ifelse(data2$pvalue < 0.05,
                       ifelse(data2$pvalue < 0.01,"**","*"),
                       "")
data2$pvalue[1:20]
data2$type <- data2$cor


summary(data2)
data3 <- data2[order(data2$immuneGene,data2$cor),]
head(data3)
dim(data3)
data4 <- data3[data3$pvalue < 0.05,]
dim(data4)
summary(data4)

p <- ggplot(data4,aes(x=gene,y=immuneGene)) +
  geom_point(aes(colour = cor, size=pvalue)) +
  labs(x="",y="Cholestasis genes")
p <- p + scale_colour_gradient2(low = "blue", high = "red", mid = "white",
                                midpoint = 0, limit = c(-1, 1), space = "Lab",
                                name="Pearson\nCorrelation")
p <- p + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  theme(axis.text=element_text(size = 15)) +
  theme(axis.text.x=element_text(colour = "black",angle=0,hjust=0.5,size = 15)) +
  theme(axis.text.y=element_text(colour = "black", vjust=0,size = 15)) +
  theme(axis.title =element_text(size = 20)) +
  theme(text = element_text(size = 15))

p+rotate_x_text(45)
ggsave("Cholestasis.pdf")


library(ggplot2)
library(ggExtra)
rt=read.table("symbol.txt",sep="\t",header=T,check.names=F,row.names = 1)
rt=log2(rt+1)
dat<-as.data.frame(t(rt))


corr_eqn <- function(x,y,digits=3) {
  test <- cor.test(x,y,type="pearson")
  paste(paste0("n = ",length(x)),
        paste0("r = ",round(test$estimate,digits),"(Pearson)"),
        paste0("p.value= ",round(test$p.value,digits)),
        sep = ", ")
}
gene<-as.numeric(dat$IL32)
imucell<-dat$ABCC2
corr_eqn(gene,imucell)


gg<-ggplot(dat, aes(x=gene, y=imucell)) + 
  geom_point(color = "black") + 
  geom_smooth(method="loess", se=F,color="blue") + 
  
  
  labs( 
    y="ABCC2", 
    x="IL32", 
    title="Scatterplot")+ 
  
  labs(title = paste0(corr_eqn(gene,imucell)))+
  theme_bw()
gg
gg2 <- ggMarginal(gg, type="density")
gg2 <- ggMarginal(gg, type="density",xparams = list(fill ="orange"),
                  yparams = list(fill ="skyblue"))
