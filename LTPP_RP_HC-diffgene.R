############# 2020-11-20 ####################
library(edgeR)
library("limma")
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(ggsignif)
library("ggpubr")
################################ 入组的83个样本 ########################
df<-read.table("featureCounts-CY_ZY_HC-wangj.txt",header=T,sep="\t",row.name=1)
names(df)
raw_count<-df[,-15]   ### CY-ZY 配对 12 对 ，HC 13 个样本
names(raw_count)   
group_list <- c(rep("LTPP",12),rep("RP",12),rep("HC",13))
##筛选前
raw_count_before_filter <- DGEList(raw_count) #构建 DGEList类变量
raw_count_before_filter$samples$group <- group_list

##根据表达量筛选 cpm >= 0.5
ncpm_before_filter<-cpm(raw_count_before_filter)
expr<-(rowSums(ncpm_before_filter[,1:12] >= 0.5)>=6 | rowSums(ncpm_before_filter[,13:24] >=0.5 )>=6 | rowSums(ncpm_before_filter[,25:37] >=0.5)>=7)
raw_count_after_filter<-raw_count_before_filter[expr,]
dim(raw_count_after_filter)   # 21056
#筛选后 标准化 edgeR 在默认情况下，执行TMM标准化程序以考虑样本之间的不同测序深度
raw_count_after_DEGlist <- calcNormFactors(raw_count_after_filter, method = "TMM")
ncpm_after_filter <- cpm(raw_count_after_DEGlist)
lcpm<-cpm(raw_count_after_DEGlist,log=TRUE,prior.count=2)
write.table(ncpm_after_filter, "1-total-ncmp-0.5.txt")   

##(2)获取转录基因的信息（基于gtf文件）
library(dplyr)
gtf1 <- rtracklayer::import('../Homo_sapiens.GRCh38.84.gtf')
gtf_df <- as.data.frame(gtf1)
test <- gtf_df[1:20,]
View(test)
#选取gene_name,gene_id,gene_biotype,start,end,width,strand
geneid_df <- dplyr::select(gtf_df,c(gene_name,gene_id,type,gene_biotype,start,end,width,strand))
sort( table( geneid_df$gene_biotype ) )
##去重复
index <- duplicated(geneid_df[,2])
geneid_df <- geneid_df[!index, ]

##与基因表达的ncpm值合并
gene_cpm<-read.table("1-total-ncmp-0.5.txt",header=T)
gene_inf<-merge(geneid_df,gene_cpm,by="gene_id")
write.table(gene_inf,"2-total_gene_information.txt",row.names=F)

### 每组样本的平均基因表达量 ##############################

LTPP_mean<-gene_inf[,9:20] %>% mutate(means=rowMeans(.))
RP_mean<-gene_inf[,21:32] %>% mutate(means=rowMeans(.))
HC_mean<-gene_inf[,33:45] %>% mutate(means=rowMeans(.))
total<-cbind(LTPP_mean$means,RP_mean$means)
total<-cbind(total,HC_mean$means)
rownames(total)<-gene_inf$gene_name
colnames(total)<-c("LTPP","RP","HC")
write.table(total,"1-gene-total-mean.txt")

##### PCA ##################
library(factoextra)
data<-t(lcpm)##转换数据至行为sample,列为gene
data<-as.data.frame(data)##注意数据要转换为数据框
res.pca <- prcomp(data, scale = TRUE)
opar<-par(no.readonly=TRUE)
pdf(".\\figure\\0-PCA-0.5.pdf",width=8,height=8)
fviz_pca_ind(res.pca,title="",
             col.ind = group_list, # 颜色对应group信息
             palette = c(brewer.pal(11,"RdGy")[9],brewer.pal(8,"Reds")[8],brewer.pal(8,"YlOrBr")[5]),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             #ellipse.type = "convex",
             legend.title = "Groups",label="none",## Legend名称
             repel = TRUE,pointsize = 3,labelsize = 6, ggtheme=theme_bw(base_size = 44),font.legend = c(25, "plain", "black"),
             font.main = c(25, "plain", "black"),font.y = c(25, "plain", "black"),font.x = c(25, "plain", "black"),  
             font.tickslab=c(25, "plain", "black")
)
dev.off()

###对样本进行无监督聚类
display.brewer.all()
lcpm <- cpm(raw_count_after_DEGlist, log=TRUE)
col.group <- factor(group_list,levels=c("LTPP","RP","HC"))
brewer.pal(8,"Reds")[8],brewer.pal(8,"Blues")[7],brewer.pal(8,"Greens")[7],
levels(col.group) <-c(brewer.pal(8,"Reds")[8],brewer.pal(8,"YlOrBr")[5],brewer.pal(11,"RdGy")[9])
col.group <- as.character(col.group)
opar<-par(no.readonly=TRUE)
pdf(".\\figure\\1-raw_count_plotMDS-无监督聚类.pdf",width=9,height=8)
par(mai=c(1,1,1,1))
plotMDS(lcpm,labels=colnames(lcpm),col=col.group,cex=1.5,cex.lab=2,cex.axis=2,pch=group_list)
dev.off()

##(3)总表达量、以及biotype前五箱线图
##图一： 总表达量
gene_inf<-read.table("2-total_gene_information.txt",header=T)
ncpm<-gene_inf[,c(1,9:45)]
head(ncpm)
data_m <- melt(ncpm)
data_m$group<-c(rep("LTPP",nrow(ncpm)*12),rep("RP",nrow(ncpm)*12),rep("HC",nrow(ncpm)*13))
data_m$group<-factor(data_m$group,levels=c("LTPP","RP","HC"))
col.group <-c(LTPP=brewer.pal(8,"Reds")[8],RP=brewer.pal(8,"YlOrBr")[5],HC=brewer.pal(11,"RdGy")[9])
options(scipen = 2000,digits=0)
pdf(".\\figure\\2-total-cpm-0.5.pdf",width=7,height=8)
ggplot(data_m,aes(x=group,y=value,fill=group))+
     geom_boxplot(width=.8,size=1)+ 
     #stat_compare_means(method = "t.test",label = "p.signif",size=8,ref.group = "HC")+
     scale_fill_manual(values=col.group)+   
     scale_y_log10()+
     theme_bw()+
     ylab('CPM')+#修改y轴名称
     xlab('')+#修改x轴名称
     theme(title=element_text(size=20),axis.text.x = element_text(size=20,angle = 45,hjust = 1),axis.text.y = element_text(size=20),
     axis.text=element_text(color="black"),panel.background=element_rect(color="black",size=2),plot.margin=unit(rep(3,4),'lines'),
     legend.text = element_text(size = 20),legend.key.size=unit(1,'cm'),panel.grid=element_blank())
dev.off()

##图二：前五种gene_type
biotype_5<-sort(table(gene_inf$gene_biotype),decreasing=TRUE)[1:5]
gene_biotype<-gene_inf[which(gene_inf$gene_biotype == 'protein_coding' | gene_inf$gene_biotype == 'antisense' | gene_inf$gene_biotype == 'lincRNA'| gene_inf$gene_biotype == 'processed_pseudogene' | gene_inf$gene_biotype == 'sense_intronic'),]
names(gene_biotype)
gene_biotype<-gene_biotype[,c(1,4,9:45)]
head(gene_biotype)
biotype_m<-melt(gene_biotype)
biotype_m$group<-c(rep("LTRR",nrow(gene_biotype)*12),rep("RP",nrow(gene_biotype)*12),rep("HC",nrow(gene_biotype)*13))
biotype_m$group<-factor(biotype_m$group,levels=c("LTRR","RP","HC"))
biotype_m$gene_biotype<-factor(biotype_m$gene_biotype,levels=c("protein_coding","antisense","lincRNA","processed_pseudogene","sense_intronic"))
names(biotype_m)[4]<-"CPM"
options(scipen = 2000,digits=0)
col.group <-c(LTRR=brewer.pal(8,"Reds")[8],RP=brewer.pal(8,"YlOrBr")[5],HC=brewer.pal(11,"RdGy")[9])
pdf(".\\figure\\3-total-biotype.pdf",width=22,height=16)
ggplot(biotype_m,aes(x=gene_biotype,y=CPM,color=group),palette = "jco")+
     geom_boxplot(size=1.5)+ 
     stat_compare_means(method = "anova",label = "p.signif",size=16)+  
     scale_y_log10()+
     scale_color_manual(values=col.group)+
     theme_bw()+
     ylab('CPM')+#修改y轴名称
     xlab('')+#修改x轴名称  
     labs(color='Gene type')+  
     scale_x_discrete(breaks = c("protein_coding","antisense","lincRNA","processed_pseudogene","sense_intronic"),
                      label = c("Protein\ncoding","Antisense","LincRNA","Processed\npseudogene","Sense\nintronic"))+
     theme(title=element_text(size=40),axis.text.x = element_text(size=40),axis.text.y = element_text(size=40),
     axis.text=element_text(color="black"),panel.background=element_rect(color="black",size=4),plot.margin=unit(rep(1,4),'lines'),
     legend.text = element_text(size = 30),legend.key.size=unit(2,'cm'),panel.grid=element_blank())
dev.off()

## LTRR、RP、HC不同类型转录基因的个数
LTPP_gene<-gene_inf[,c(1,2,4,9:20)]
LTPP_gene$gene_biotype<-factor(LTPP_gene$gene_biotype)
LTPP_noexpr<-rowSums(LTPP_gene[,4:ncol(LTPP_gene)] >= 0.5 )>=6
LTPP_gene_filted <- LTPP_gene[LTPP_noexpr,]
dim(LTPP_gene_filted)
summary(LTPP_gene_filted$gene_biotype)

RP_gene<-gene_inf[,c(1,2,4,21:32)]
RP_gene$gene_biotype<-factor(RP_gene$gene_biotype)
RP_noexpr<-rowSums(RP_gene[,4:ncol(RP_gene)] >=0.5 )>=6
RP_gene_filted <- RP_gene[RP_noexpr,]
dim(RP_gene_filted)
summary(RP_gene_filted$gene_biotype)

HC_gene<-gene_inf[,c(1,2,4,33:45)]
HC_gene$gene_biotype<-factor(HC_gene$gene_biotype)
HC_noexpr<-rowSums(HC_gene[,4:ncol(HC_gene)] >= 0.5 )>=7
HC_gene_filted <- HC_gene[HC_noexpr,]
dim(HC_gene_filted)
summary(HC_gene_filted$gene_biotype)

library(VennDiagram)
venn.diagram(list(RP=RP_gene_filted$gene_id,LTPP=LTPP_gene_filted$gene_id,HC=HC_gene_filted$gene_id),
             resolution = 500, alpha=c(0.5,0.5,0.5),rotation.degree=180,,cat.pos=c(120, 240, 0),
             fill=c(brewer.pal(8,"YlOrBr")[5],brewer.pal(8,"Reds")[8],brewer.pal(11,"RdGy")[9]), 
             cex = 2,cat.cex = 2 ,cat.dist = 0.05,lty = "blank",margin = 0.05,
             filename = ".\\figure\\VennDiagram_specific.png"
)


### （4）差异基因分析
#1.创建设计矩阵和对比
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(raw_count_after_DEGlist)
cont.matrix  = makeContrasts(contrasts=c('LTPP-RP','LTPP-HC','RP-HC'),levels=design)

#2.找差异基因
v <- voom(raw_count_after_DEGlist,design,plot=TRUE)
vfit <- lmFit(v,design)
vfit=contrasts.fit(vfit,cont.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
#3.检查DE基因的数量
summary(decideTests(efit,lfc=1,p.value =0.01))
summary(decideTests(efit,lfc=2,p.value =0.05))
summary(decideTests(efit,lfc=1.5,p.value =0.05))

#4.导出所有的有显著差异的基因
### LTPP_RP
tmpOutput1=topTable(efit,coef='LTPP-RP',n=Inf)
all_voom1 = na.omit(tmpOutput1)
all_diff1 <- all_voom1[abs(all_voom1$logFC)>=1.5&all_voom1$adj.P.Val<0.05,]  ##条件需要根据后续要求来改变
dim(all_diff1)
write.table(all_voom1,"3-total-lr-default-gene.txt")
write.table(all_diff1,"3-total-lr-diffgene_FC.txt")

### LTPP_HC
tmpOutput2=topTable(efit,coef='LTPP-HC',n=Inf)
all_voom2 = na.omit(tmpOutput2)
all_diff2 <- all_voom2[abs(all_voom2$logFC)>=1.5&all_voom2$adj.P.Val<0.05,]  ##条件需要根据后续要求来改变
dim(all_diff2)
write.table(all_voom2,"3-total-lh-default-gene.txt")
write.table(all_diff2,"3-total-lh-diffgene_FC.txt")

### RP_HC
tmpOutput3=topTable(efit,coef='RP-HC',n=Inf)
all_voom3 = na.omit(tmpOutput3)
all_diff3 <- all_voom3[abs(all_voom3$logFC)>=1.5&all_voom3$adj.P.Val<0.05,]  ##条件需要根据后续要求来改变
dim(all_diff3)
write.table(all_voom3,"3-total-rh-default-gene.txt")
write.table(all_diff3,"3-total-rh-diffgene_FC.txt")

## 提取三组差异基因的基因信息
gene_inf<-read.table("2-total_gene_information.txt",header=T)
all_diff1<-read.table("3-total-lr-diffgene_FC.txt",header=T)
all_default1<-read.table("3-total-lr-default-gene.txt",header=T)
all_diff1$gene_id<-rownames(all_diff1)
all_default1$gene_id<-rownames(all_default1)
LTPP_RP_gene_inf<-merge(x=all_diff1,y=gene_inf,by="gene_id",all.x = TRUE)
LTPP_RP_default_inf<-merge(x=all_default1,y=gene_inf,by="gene_id",all.x = TRUE)
dim(LTPP_RP_gene_inf)
write.table(LTPP_RP_gene_inf,"4-total-lr_gene_inf.txt",row.names=F)
write.table(LTPP_RP_default_inf,"4-total-lr_default_inf.txt",row.names=F)

all_diff2<-read.table("3-total-lh-diffgene_FC.txt",header=T)
all_default2<-read.table("3-total-lh-default-gene.txt",header=T)
all_diff2$gene_id<-rownames(all_diff2)
all_default2$gene_id<-rownames(all_default2)
LTPP_HC_gene_inf<-merge(x=all_diff2,y=gene_inf,by="gene_id",all.x = TRUE)
LTPP_HC_default_inf<-merge(x=all_default2,y=gene_inf,by="gene_id",all.x = TRUE)
dim(LTPP_HC_gene_inf)
write.table(LTPP_HC_gene_inf,"4-total-lh_gene_inf.txt",row.names=F)
write.table(LTPP_HC_default_inf,"4-total-lh_default_inf.txt",row.names=F)

all_diff3<-read.table("3-total-rh-diffgene_FC.txt",header=T)
all_default3<-read.table("3-total-rh-default-gene.txt",header=T)
all_diff3$gene_id<-rownames(all_diff3)
all_default3$gene_id<-rownames(all_default3)
RP_HC_gene_inf<-merge(x=all_diff3,y=gene_inf,by="gene_id",all.x = TRUE)
RP_HC_default_inf<-merge(x=all_default3,y=gene_inf,by="gene_id",all.x = TRUE)
dim(RP_HC_gene_inf)
write.table(RP_HC_gene_inf,"4-total-rh_gene_inf.txt",row.names=F)
write.table(RP_HC_default_inf,"4-total-rh_default_inf.txt",row.names=F)


##图四 火山图
library(ggplot2)
library(ggrepel)
#读取数据
Dat<-read.table('4-total-lr_default_inf.txt',header = T)
Dat<-Dat[,c(8,2:7)]
#作图
#确定是上调还是下调，用于给图中点上色）
Dat$threshold = factor(ifelse(Dat$adj.P.Val < 0.05 & abs(Dat$logFC) >= 1.5, ifelse(Dat$logFC>= 1.5 ,'Up','Down'),'NoSign'),levels=c('Up','Down','NoSign'))
sign_data<-Dat[Dat$adj.P.Val<0.05 & abs(Dat$logFC)>1.5,]
down_10 <- sign_data[order(sign_data$logFC),][1:10,]
up_10 <- sign_data[order(sign_data$logFC,decreasing=TRUE),][1:10,]
my_geom_text<-rbind(down_10,up_10)
pdf(".\\figure\\6-volcano_LTPP_RP.pdf",width=6,height=7.5)
ggplot(Dat,aes(x=logFC,y=-log10(adj.P.Val),color=threshold))+
  geom_point()+
  scale_color_manual(values=c(brewer.pal(8,"Reds")[8],brewer.pal(8,"Blues")[8],"#808080"))+ #确定点的颜色
  theme_bw()+#修改图片背景
  guides(color=FALSE)+ #不显示图例
  theme(axis.text=element_text(size=20,color="black"),
        axis.title=element_text(size=25,color="black"),
        title= element_text(size=28,color="black"),
        #panel.grid=element_blank(),
        plot.title = element_text(hjust = 0.5), 
        plot.margin=unit(c(1,1,1,1),'lines'))+ 
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  labs(title="LTPP-RP")+
  geom_vline(xintercept=c(-1.5,1.5),lty=3,col="black",lwd=1) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=1) + #添加竖线padj<0.05
  geom_label_repel(
    data = my_geom_text,
    aes(label = gene_name),box.padding=unit(0.4, "lines"),force = 10,
    size = 5,color="black",#fontface="bold",
    segment.color = "black", show.legend = FALSE )#添加关注的点的基因名
dev.off()

Dat<-read.table('4-total-lh_default_inf.txt',header = T)
Dat<-Dat[,c(8,2:7)]
Dat$threshold = factor(ifelse(Dat$adj.P.Val < 0.05 & abs(Dat$logFC) >= 1.5, ifelse(Dat$logFC>= 1.5 ,'Up','Down'),'NoSign'),levels=c('Up','Down','NoSign'))
sign_data<-Dat[Dat$adj.P.Val<0.05 & abs(Dat$logFC)>1.5,]
down_10 <- sign_data[order(sign_data$logFC),][1:10,]
up_10 <- sign_data[order(sign_data$logFC,decreasing=TRUE),][1:10,]
my_geom_text<-rbind(down_10,up_10)

pdf(".\\figure\\6-volcano_LTPP_hc.pdf",width=6,height=7.5)
ggplot(Dat,aes(x=logFC,y=-log10(adj.P.Val),color=threshold))+
  geom_point()+
  scale_color_manual(values=c(brewer.pal(8,"Reds")[8],brewer.pal(8,"Blues")[8],"#808080"))+ #确定点的颜色
  theme_bw()+#修改图片背景
  guides(color=FALSE)+ #不显示图例
  theme(axis.text=element_text(size=20,color="black"),
        axis.title=element_text(size=25,color="black"),
        title= element_text(size=28,color="black"),
        plot.title = element_text(hjust = 0.5),plot.margin=unit(c(1,1,1,1),'lines'))+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  labs(title="LTPP-HC")+
  geom_vline(xintercept=c(-1.5,1.5),lty=3,col="black",lwd=1) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=1) + #添加竖线padj<0.05
  geom_label_repel(
    data = my_geom_text,
    aes(label = gene_name),box.padding=unit(0.4, "lines"),force = 10,
    size = 5,color="black",#fontface="bold",
    segment.color = "black", show.legend = FALSE )#添加关注的点的基因
dev.off()

Dat<-read.table('4-total-rh_default_inf.txt',header = T)
Dat<-Dat[,c(8,2:7)]
Dat$threshold = factor(ifelse(Dat$adj.P.Val < 0.05 & abs(Dat$logFC) >= 1.5, ifelse(Dat$logFC>= 1.5 ,'Up','Down'),'NoSign'),levels=c('Up','Down','NoSign'))
sign_data<-Dat[Dat$adj.P.Val<0.05 & abs(Dat$logFC)>1.5,]
down_10 <- sign_data[order(sign_data$logFC),][1:10,]
up_10 <- sign_data[order(sign_data$logFC,decreasing=TRUE),][1:10,]
my_geom_text<-rbind(down_10,up_10)

pdf(".\\figure\\6-volcano_RP_hc.pdf",width=6,height=7.5)
ggplot(Dat,aes(x=logFC,y=-log10(adj.P.Val),color=threshold))+
  geom_point()+
  scale_color_manual(values=c(brewer.pal(8,"Reds")[8],brewer.pal(8,"Blues")[8],"#808080"))+ #确定点的颜色
  theme_bw()+#修改图片背景
  guides(color=FALSE)+ #不显示图例
  theme(axis.text=element_text(size=20,color="black"),
        axis.title=element_text(size=25,color="black"),
        title= element_text(size=28,color="black"),
        plot.title = element_text(hjust = 0.5),plot.margin=unit(c(1,1,1,1),'lines'))+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  labs(title="RP-HC")+
  geom_vline(xintercept=c(-1.5,1.5),lty=3,col="black",lwd=1) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=1) + #添加竖线padj<0.05
  geom_label_repel(
    data = my_geom_text,
    aes(label = gene_name),box.padding=unit(0.4, "lines"),force = 10,
    size = 5,color="black",#fontface="bold",
    segment.color = "black", show.legend = FALSE )#添加关注的点的基因名
dev.off()


## 提取所有3个组的差异基因
LTPP_RP_diffname<-rownames(all_diff1)
LTPP_HC_diffname<-rownames(all_diff2)
RP_HC_diffname<-rownames(all_diff3)
diffname<-c(LTPP_RP_diffname,LTPP_HC_diffname,RP_HC_diffname)
unique_diffname<-unique(diffname)
diffgene<-gene_inf[gene_inf$gene_id %in% unique_diffname,]
write.table(diffgene,"5-diff_gene-total-inf.txt",row.names=F)
names(diffgene)
dim(diffgene)
diff_gene_ncpm<-diffgene[,c(1,9:45)]

library(pheatmap)
annotation_col <- data.frame(SampleType = factor(c(rep("LTPP",12),rep("RP",12),rep("HC",13))))  #增加Time，CellType分组信息
rownames(annotation_col) = colnames(diff_gene_ncpm)[-1]   
ann_colors = list(SampleType = c(LTPP=brewer.pal(8,"Reds")[8],RP=brewer.pal(8,"YlOrBr")[5],HC=brewer.pal(11,"RdGy")[9]))
#bk = unique(c(seq(0,3, length=100)))
bk = unique(c(seq(-3,3, length=100)))
pdf(".\\figure\\5-LTPP_RP_HC-pheatmap.pdf",width=12,height=16) 
p=pheatmap(log10(diff_gene_ncpm[,-1]+1),scale = "row",show_rownames = FALSE,show_colnames = F,clustering_method="average",
         cellwidth=12,cellheight=0.1,treeheight_row=0,cluster_col =F,fontsize=15,annotation_col = annotation_col,annotation_colors = ann_colors,
         color = colorRampPalette(colors = c(brewer.pal(8,"Blues")[8],"white",brewer.pal(8,"Reds")[8]))(100),
         breaks=bk,annotation_names_col=FALSE,cutree_rows=3
        )
dev.off()
row_cluster1 <- cutree(p$tree_row,k=3)
newOrder1 <- diff_gene_ncpm[p$tree_row$order,]
cluster<-row_cluster1[match(rownames(newOrder1),names(row_cluster1))]
newOrder1<-cbind(newOrder1,cluster)
newOrder1<-as.data.frame(newOrder1)
newOrder1$cluster<-factor(newOrder1$cluster)
summary(newOrder1$cluster)          ### 从上至下：cluster3=1891(RP),cluster2=1797(HC),cluster1=2152

LTPP_mean<-diff_gene_ncpm[,2:13] %>% mutate(means=rowMeans(.))
RP_mean<-diff_gene_ncpm[,14:25] %>% mutate(means=rowMeans(.))
HC_mean<-diff_gene_ncpm[,26:38] %>% mutate(means=rowMeans(.))
total<-cbind(LTPP_mean$means,RP_mean$means)
total<-cbind(total,HC_mean$means)
rownames(total)<-diff_gene_ncpm$gene_id
colnames(total)<-c("LTPP","RP","HC")
write.table(total,"5-diff_gene-total-mean.txt")
bk = unique(c(seq(-1,1, length=100)))
p=pheatmap(log10(total+1),scale = "row",show_rownames = FALSE,cluster_col =F,fontsize=18,clustering_distance_rows="canberra",
        color = colorRampPalette(colors = c(brewer.pal(8,"Blues")[8],"white",brewer.pal(8,"Reds")[8]))(100),breaks=bk,
        legend_breaks=c(-1,1),angle_col=45,cellwidth=50,cellheight=0.06,cutree_rows=3,
        )
row_cluster <- cutree(p$tree_row,k=3)
## 对结果进行排序
newOrder <- total[p$tree_row$order,]
cluster<-row_cluster[match(rownames(newOrder),names(row_cluster))]
newOrder<-cbind(newOrder,cluster)
newOrder<-as.data.frame(newOrder)
newOrder$cluster<-factor(newOrder$cluster)
summary(newOrder$cluster)          ### 从上至下：cluster3=1891(RP),cluster2=1797(HC),cluster1=2152
write.table(newOrder,"5-pheatmap-cluster.txt")
### 添加基因名字 #############
df<-read.table("5-pheatmap-cluster.txt",header=T)
gene_inf<-read.table("2-total_gene_information.txt",header=T)
df_name<-merge(df,gene_inf[,1:2],by="gene_id",all.x=T)
write.table(df_name,"5-pheatmap-cluster-genename.txt")
### 
pheatmap_new<-read.table("5-pheatmap-cluster.txt",header=T,row.names=1)
annotation_row <- data.frame(SampleType = factor(c(rep("LTPP",2152),rep("RP",1891),rep("HC",1797))))  #增加Time，CellType分组信息
rownames(annotation_row) = rownames(pheatmap_new)  
ann_colors = list(SampleType = c(LTPP=brewer.pal(8,"Reds")[8],RP=brewer.pal(8,"YlOrBr")[5],HC=brewer.pal(11,"RdGy")[9]))
pdf(".\\figure\\5-LTPP_RP_HC-pheatmap-mean.pdf",width=6,height=8) 
pheatmap(log10(pheatmap_new[,-4]+1),scale = "row",show_rownames = FALSE,cluster_col =F,fontsize=18,cluster_row =T,
        clustering_distance_rows="canberra",treeheight_row=0, #clustering_method="ward.D2",#"average",#"complete",#,"correlation""maximum""canberra""binary""euclidean""ward.D2"
        color = colorRampPalette(colors = c(brewer.pal(8,"Blues")[8],"white",brewer.pal(8,"Reds")[8]))(100),breaks=unique(c(seq(-1,1, length=100))),
        annotation_row = annotation_row,annotation_colors = ann_colors,annotation_names_row=FALSE,
        legend_breaks=c(-1,1),angle_col=45,cellwidth=50,cellheight=0.06 ,cutree_rows=3,
        legend_position="left"
        )
dev.off()

######### specific ##################
###筛选每组特异表达的基因
gene_inf<-read.table("2-total_gene_information.txt",header=T)
LTPP_RP_default_inf<-read.table("4-total-lr_default_inf.txt",header=T)
LTPP_HC_default_inf<-read.table("4-total-lh_default_inf.txt",header=T)
Mild_HC_default_inf<-read.table("4-total-rh_default_inf.txt",header=T)
dim(LTPP_RP_default_inf)
names(gene_inf)

##长阳组
LTPP_onlyexpr<-(rowSums(gene_inf[,9:20] >= 0.5 )>=6 & rowSums(gene_inf[,21:32] >= 0.5 )<6 & rowSums(gene_inf[,33:45] >= 0.5 )<7 )
LTPP_only <- gene_inf[LTPP_onlyexpr,]
dim(LTPP_only)
LTPP_only$gene_biotype<-factor(LTPP_only$gene_biotype)
summary(LTPP_only$gene_biotype)
LTPP_only_inf<-merge(y=LTPP_only,x=LTPP_RP_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(LTPP_only_inf)[2:3]<-c("LR_logFC","LR_adj.p.value")
LTPP_only_inf<-merge(y=LTPP_only_inf,x=LTPP_HC_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(LTPP_only_inf)[2:3]<-c("S=LH_logFC","LH_adj.p.value")
LTPP_only_inf<-merge(y=LTPP_only_inf,x=RP_HC_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(LTPP_only_inf)[2:3]<-c("RH_logFC","RH_adj.p.value")
write.table(LTPP_only_inf,"6-LTPP-specific.txt",row.names=F)
dim(LTPP_only_inf)

##轻症组
RP_onlyexpr<-(rowSums(gene_inf[,9:20] >= 0.5 )<6 & rowSums(gene_inf[,21:32] >= 0.5 )>=6 & rowSums(gene_inf[,33:45] >= 0.5 )<7 )
RP_only <- gene_inf[RP_onlyexpr,]
dim(RP_only)
RP_only$gene_biotype<-factor(RP_only$gene_biotype)
summary(RP_only$gene_biotype)
RP_only_inf<-merge(y=RP_only,x=LTPP_RP_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(RP_only_inf)[2:3]<-c("LR_logFC","LR_adj.p.value")
RP_only_inf<-merge(y=RP_only_inf,x=LTPP_HC_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(RP_only_inf)[2:3]<-c("LH_logFC","LH_adj.p.value")
RP_only_inf<-merge(y=RP_only_inf,x=RP_HC_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(RP_only_inf)[2:3]<-c("RH_logFC","RH_adj.p.value")
write.table(RP_only_inf,"6-RP-specific.txt",row.names=F)
dim(RP_only_inf)

##健康组
HC_onlyexpr<-(rowSums(gene_inf[,9:20] >= 0.5 )<6 & rowSums(gene_inf[,21:32] >= 0.5 )<6 & rowSums(gene_inf[,33:45] >= 0.5 )>=7 )
HC_only <- gene_inf[HC_onlyexpr,]
dim(HC_only)
HC_only$gene_biotype<-factor(HC_only$gene_biotype)
summary(HC_only$gene_biotype)
HC_only_inf<-merge(y=HC_only,x=LTPP_RP_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(HC_only_inf)[2:3]<-c("LR_logFC","LR_adj.p.value")
HC_only_inf<-merge(y=HC_only_inf,x=LTPP_HC_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(HC_only_inf)[2:3]<-c("LH_logFC","LH_adj.p.value")
HC_only_inf<-merge(y=HC_only_inf,x=RP_HC_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(HC_only_inf)[2:3]<-c("RH_logFC","RH_adj.p.value")
write.table(HC_only_inf,"6-HC-specific.txt",row.names=F)
dim(HC_only_inf)

######### pattern ##################
diffgene<-read.table("5-diff_gene-total-inf.txt",header=T)
diffgene1<-merge(y=diffgene,x=LTPP_RP_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(diffgene1)[c(2,3)]<-c("LR_logFC","LR_adj.p.value")
diffgene2<-merge(y=diffgene1,x=LTPP_HC_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(diffgene2)[2:3]<-c("LH_logFC","LH_adj.p.value")
diffgene3<-merge(y=diffgene2,x=RP_HC_default_inf[,c(1,2,6)],by="gene_id",all.y = TRUE)
names(diffgene3)[2:3]<-c("RH_logFC","RH_adj.p.value")
dim(diffgene3)
write.table(diffgene3,"7-diff_gene_patternFC.txt",row.names=F)

########################## HML 热图 #############################
library(pheatmap)
HML<-read.table("7-HML-top30.txt",header=T,sep="\t",row.names=1)
pdf(".\\figure\\7-HML-top30.pdf",width=5,height=10)
pheatmap(log10(HML+1),cellwidth=20,cellheight=20,treeheight_row=0,cluster_col =F,cluster_row =F,scale="row",
         fontsize=16,legend_breaks=c(-1,1),breaks=unique(c(seq(-1,1,length=90))),angle_col=45,
         color = colorRampPalette(colors = c(brewer.pal(8,"BrBG")[7],"white",brewer.pal(8,"RdBu")[1]))(90)
)
dev.off()

LMH<-read.table("7-LMH-top30.txt",header=T,sep="\t",row.names=1)
pdf(".\\figure\\7-LMH-top30.pdf",width=5,height=10)
pheatmap(log10(LMH+1),cellwidth=20,cellheight=20,treeheight_row=0,cluster_col =F,cluster_row =F,scale="row",
         fontsize=16,legend_breaks=c(-1,1),breaks=unique(c(seq(-1,1,length=90))),angle_col=45,
         color = colorRampPalette(colors = c(brewer.pal(8,"BrBG")[7],"white",brewer.pal(8,"RdBu")[1]))(90)
)
dev.off()

######## 小提琴图  ###########################
##小提琴图
HML1<-read.table("7-HML-top30.txt",header=T,sep="\t")
HML_violin<-melt(HML1)
names(HML_violin)<-c("gene","group","cpm")
HML_violin$group<-factor(HML_violin$group,levels=c("LTPP","RP","HC"))
compaired <- list(c("LTPP","RP"),c("LTPP","HC"),c("RP","HC"))
pdf(".\\figure\\7-HML_point.pdf",width=4,height=6)
ggplot(HML_violin,aes(x=group,y=log10(cpm+1),color=group,fill=group))+
     #geom_violin(width=1,alpha=0.5)+
     geom_point(position=position_jitter(width=.2, height=0),alpha=0.5,size=4)+
     geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T,test = t.test,textsize = 9,color="black",size=1)+
     scale_color_manual(values=c(brewer.pal(8,"Reds")[8],brewer.pal(8,"YlOrBr")[5],brewer.pal(11,"RdGy")[9]))+  
     scale_fill_manual(values=c(brewer.pal(8,"Reds")[8],brewer.pal(8,"YlOrBr")[5],brewer.pal(11,"RdGy")[9]))+
     scale_y_continuous(limits=c(0, 6))+
     theme_bw()+guides(fill=FALSE,color=FALSE)+ #不显示图例
     ylab('logCPM')+#修改y轴名称
     xlab('')+#修改x轴名称
     theme(title=element_text(size=25,face= "bold"),axis.text.x = element_text(size=25,angle = 45,hjust = 1,face= "bold"),axis.text.y = element_text(size=25,face= "bold"),
     axis.text=element_text(color="black",face= "bold"),plot.margin=unit(rep(1,4),'lines'),
     panel.border = element_blank(),axis.line = element_line(colour = "black",size=1.5), ###去掉边框，并把坐标轴显示加粗
     panel.grid=element_blank())
dev.off()

LMH1<-read.table("7-LMH-top30.txt",header=T,sep="\t")
LMH_violin<-melt(LMH1)
names(LMH_violin)<-c("gene","group","cpm")
LMH_violin$group<-factor(LMH_violin$group,levels=c("LTPP","RP","HC"))
compaired <- list(c("LTPP","RP"),c("LTPP","HC"),c("RP","HC"))
pdf(".\\figure\\7-LMH_point.pdf",width=4,height=6)
ggplot(LMH_violin,aes(x=group,y=log10(cpm+1),color=group,fill=group))+
     #geom_violin(width=2,alpha=0.5)+
     geom_point(position=position_jitter(width=.2, height=0),alpha=0.5,size=4)+
     geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T,test = t.test,textsize = 9,color="black",size=1)+
     scale_color_manual(values=c(brewer.pal(8,"Reds")[8],brewer.pal(8,"YlOrBr")[5],brewer.pal(11,"RdGy")[9]))+  
     scale_fill_manual(values=c(brewer.pal(8,"Reds")[8],brewer.pal(8,"YlOrBr")[5],brewer.pal(11,"RdGy")[9]))+
     scale_y_continuous(limits=c(0, 4))+
     theme_bw()+guides(fill=FALSE,color=FALSE)+ #不显示图例
     ylab('logCPM')+#修改y轴名称
     xlab('')+#修改x轴名称
     theme(title=element_text(size=25,face= "bold"),axis.text.x = element_text(size=25,angle = 45,hjust = 1,face= "bold"),axis.text.y = element_text(size=25,face= "bold"),
     axis.text=element_text(color="black",face= "bold"),plot.margin=unit(rep(1,4),'lines'),
     panel.border = element_blank(),axis.line = element_line(colour = "black",size=1.5), ###去掉边框，并把坐标轴显示加粗
     panel.grid=element_blank())
dev.off()

############# pattern 线图 #############################
pattern1<-read.table("7-HML.txt",header=T)
dim(pattern1)
pattern1_log<-log10(pattern1[,2:4]+1)
pattern1_log<-cbind(pattern1[,1],pattern1_log)
names(pattern1_log)[1]<-"gene_name"
pattern1_melt<-melt(pattern1_log)
pattern1_melt_sort<-pattern1_melt[sort(pattern1_melt$gene_name,index.return=TRUE)$ix,]
pattern1_melt_sort$gene<-as.numeric(gl((nrow(pattern1_melt_sort)/3), 3))
pattern1_melt_sort$sample<-rep(1:3,(nrow(pattern1_melt_sort)/3))
##画图,c("HC","CY","ZY")
ngene <- max(pattern1_melt_sort$gene)
xrange <- range(pattern1_melt_sort$sample) 
yrange <- range(pattern1_melt_sort$value)
colors <- rainbow(ngene) 
pdf(".//figure//7-HML-line.pdf",width=8,height=8)
par(pin=c(4.5,4.5))
plot(xrange, yrange, type="n", xlab="HML",
  ylab="log10(CPM+1)",main="N=1083",xaxt="n",cex.lab=2,cex.main=3,cex.axis=1.5)
axis(1,seq(1, 3, 1),c("LTPP","RP","HC"),cex.axis=2) 
for (i in 1:ngene) { 
  tree <- subset(pattern1_melt_sort[,c(4,3,5)], gene==i) 
  lines(tree$sample, tree$value, type="b", lwd=1.5,
    col=colors[i])
}
dev.off()

pattern2<-read.table("7-LMH.txt",header=T)
head(pattern2)
dim(pattern2)
pattern_log<-log2(pattern2[,2:4]+1)
pattern_log<-cbind(pattern2[,1],pattern_log)
names(pattern_log)[1]<-"gene_id"
pattern_melt<-melt(pattern_log)
pattern_melt_sort<-pattern_melt[sort(pattern_melt$gene_id,index.return=TRUE)$ix,]
pattern_melt_sort$gene<-as.numeric(gl((nrow(pattern_melt_sort)/3), 3))
pattern_melt_sort$sample<-rep(1:3,(nrow(pattern_melt_sort)/3))
##画图,c("HC","CY","ZY")
ngene <- max(pattern_melt_sort$gene)
xrange <- range(pattern_melt_sort$sample) 
yrange <- range(pattern_melt_sort$value)
colors <- rainbow(ngene) 
pdf(".//figure//7-LMH-line.pdf",width=8,height=8)
par(pin=c(4.5,4.5))
plot(xrange, yrange, type="n", xlab="LMH",
  ylab="log10(CPM+1)",main="N=1450",xaxt="n",cex.lab=2,cex.main=3,cex.axis=1.5)
axis(1,seq(1, 3, 1),c("LTPP","RP","HC"),cex.axis=2) 
for (i in 1:ngene) { 
  tree <- subset(pattern_melt_sort[,c(4,3,5)], gene==i) 
  lines(tree$sample, tree$value, type="b", lwd=1.5,
    col=colors[i]) 
} 
dev.off()

