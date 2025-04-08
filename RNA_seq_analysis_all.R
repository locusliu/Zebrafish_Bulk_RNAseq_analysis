setwd("D:/DingYi_RNAseq/batch_effect_correct")
rm(list = ls())
options(stringsAsFactors = F) # 字符串转因子，读取数据时添加

library(stringr)
# 读取csv格式用read.csv（"filenames"）；此外，还需要设置参数，如第一列为行名、分割符、设置表头等

rawcount <- read.csv("raw_gene_count_triplicate_BE_corrected.csv", header=TRUE, stringsAsFactors=TRUE, row.names = 1)
# 查看表达谱
colnames(rawcount)  # 查看列名，选取了6个样本，二分组（Dex vs untreated）
rawcount[1:3,1:6]   # 查看表达矩阵前3行前6列
dim(rawcount)       # 查看表达矩阵的行列数


group <- read.table("Triplicate_group.txt",header = T,sep = "\t",quote = "\"") # quote用来指定包围字符型变量的符号，读入后自动剔除
group
colnames(group)  # 查看列名，选取所需要的列
group <- group[match(colnames(rawcount),group$samples),c("samples","condition")] # group[]矩阵取行和列，行是match(colnames(rawcount),group$samples)，列是group的两列c("samples","sample_title")  ## match(x,y)函数找出x中每个元素在y中的位置，即要处理的样本colnames(rawcount)在group$samples列中所在的行号
group_list <- str_split_fixed(group$samples,pattern  = "_", n=2)[,2]   # str_split_fixed处理字符串向量，返回的是矩阵；str_split处理字符串向量，返回的是列表；pattern  = "_", n=2代表按照“_”分割成两列；[,2] 取分割后矩阵的第2列

# 过滤低表达基因(至少在50%样本中表达)
keep <- rowSums(rawcount>0) >=floor(0.5*ncol(rawcount)) # 逻辑值向量可以求和（TRUE=1，FALSE=0);floor函数向下取整数
table(keep)  # 返回符合过滤条件的数量
filter_count <- rawcount[keep,] #提取rawcount矩阵符合keep条件的行
filter_count[1:6,1:6]   # 查看过滤后表达矩阵前4行前4列
dim(filter_count)      # 查看过滤后表达矩阵行列数  ##过滤后有21930个基因，由于PolyA尾是mRNA和部分lnRNA的重要结构，所以此处含有少量lnRNA，基因总数稍多于数据库中的codiing gene总数

library(edgeR)
express_cpm <- cpm(filter_count)
express_cpm[1:6,1:6]
save(filter_count, express_cpm, group, file = "triplicate_Data.Rdata")

library(FactoMineR)
library(factoextra)
library(corrplot)
library(pheatmap)
library(tidyverse)

dat <- log10(express_cpm+1)
dat[1:6,1:6]
dim(dat)
sampleTree <- hclust(dist(t(dat)), method = "average")
plot(sampleTree)


dat <- as.data.frame(t(dat))
dat_pca <- PCA(dat, graph = FALSE)

group_list <- group[match(group$samples,rownames(dat)), 2]
group_list

# geom.ind: point显示点，text显示文字
# palette: 用不同颜色表示分组
# addEllipses: 是否圈起来
mythe <- theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12))

p <- fviz_pca_ind(dat_pca,
                  geom.ind = "point",  #point
                  col.ind = group_list, 
                  palette = c("cyan", "#ff4d40"),
                  labelsize = 5, pointsize = 2.5, 
                  axes.linetype = "dashed", 
                  font.family = "Arial", 
                  # Don't use default Ellipses!!!!
                  # addEllipses = TRUE,
                  #  invisible='quali',
                  ellipse.border.remove = TRUE,
                  repel = TRUE, 
                  legend.title = "Groups") +
  # ADD ggforce's ellipses
  
  
  ggforce::geom_mark_ellipse(aes(fill = group_list, 
                                 color = group_list)) +
  xlim(-150, 200) + ylim (-150, 250) +
  mythe

p

## 3.样本之间的相关性-cor----
# 选择差异变化大的基因算样本相关性
exprSet <- express_cpm
exprSet = exprSet[names(sort(apply(exprSet, 1, mad),decreasing = T)[1:800]),]
dim(exprSet)

M <- cor(exprSet,method = "spearman")
M

# 构造注释条
anno <- data.frame(group=group$condition,row.names = group$samples )

group_df = data.frame(Groups=as.factor(rep(c("A543G", "WT"), c(5,5))))

rownames(group_df) <- colnames(matrix(M))

ann_colors = list(
  Groups = c(A543G="#ff4d40", WT="cyan"))

# 保存结果
cor_samples <-pheatmap(M,display_numbers = T,
                       fontsize_number=11,
                       fontfamily = "Arial",
                       annotation_col = anno,
                       annotation_colors = ann_colors,
                       fontsize = 11,cellheight = 30,
                       cellwidth = 30,cluster_rows = T,
                       cluster_cols = T)


cor_samples

# rm(list = ls())
# options(stringsAsFactors = F)


# 加载包
library(edgeR)
library(ggplot2)

# 读取基因表达矩阵信息并查看分组信息和表达矩阵数据
lname <- load(file = "triplicate_Data.Rdata")
lname

dat <- log10(express_cpm+1)

rawcount <- read.csv("raw_gene_count_triplicate_BE_corrected.csv", header=TRUE, stringsAsFactors=TRUE, row.names = 1)

group <- read.table("Triplicate_group.txt",header = T,sep = "\t",quote = "\"") # quote用来指定包围字符型变量的符号，读入后自动剔除
group
colnames(group)  # 查看列名，选取所需要的列
group <- group[match(colnames(rawcount),group$samples),c("samples","condition")] # group[]矩阵取行和列，行是match(colnames(rawcount),group$samples)，列是group的两列c("samples","sample_title")  ## match(x,y)函数找出x中每个元素在y中的位置，即要处理的样本colnames(rawcount)在group$samples列中所在的行号
#group_list <- str_split_fixed(group$samples,pattern  = "_", n=2)[,2]   # str_split_fixed处理字符串向量，返回的是矩阵；str_split处理字符串向量，返回的是列表；pattern  = "_", n=2代表按照“_”分割成两列；[,2] 取分割后矩阵的第2列


group_list <- group[match(colnames(filter_count),group$samples), 2]
group_list

comp <- unlist(strsplit("A543G_vs_WT",split = "_vs_"))
group_list <- factor(group_list,levels = comp)
group_list
table(group_list)

design <- model.matrix(~0+group_list)
rownames(design) <- colnames(filter_count)
colnames(design) <- levels(factor(group_list))
design

# 构建edgeR的DGEList对象
DEG <- DGEList(counts=filter_count, 
               group=factor(group_list))

# 归一化基因表达分布
DEG <- calcNormFactors(DEG)

# 计算线性模型的参数
DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

# 拟合线性模型
fit <- glmFit(DEG, design)

# 进行差异分析，计算LRT值
lrt <- glmLRT(fit, contrast=c(1,-1)) 

# 提取过滤差异分析结果
# 其中LR为上面计算的LRT值，FDR值为矫正后的P值
DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG)))
head(DEG_edgeR)

fc_cutoff <- 1.5
pvalue <- 0.05

# 先把默认值都设置为normal
DEG_edgeR$regulated <- "normal"

# which返回满足两个条件的位置信息，取交集
loc_up <- intersect(which( DEG_edgeR$logFC > log2(fc_cutoff) ),
                    which( DEG_edgeR$PValue < pvalue) )

loc_down <- intersect(which(DEG_edgeR$logFC < (-log2(fc_cutoff))),
                      which(DEG_edgeR$PValue<pvalue))

# 把对应位置的regulated参数设置为up和down
DEG_edgeR$regulated[loc_up] <- "up"
DEG_edgeR$regulated[loc_down] <- "down"

table(DEG_edgeR$regulated)

library(org.Dr.eg.db)
keytypes(org.Dr.eg.db)

library(clusterProfiler)
id2symbol <- bitr(rownames(DEG_edgeR), 
                  fromType = "ENSEMBL", 
                  toType = "SYMBOL", 
                  OrgDb = org.Dr.eg.db)
head(id2symbol)

DEG_edgeR <- cbind(GeneID=rownames(DEG_edgeR),DEG_edgeR)
DEG_edgeR_symbol <- merge(id2symbol,DEG_edgeR,
                          by.x="ENSEMBL",by.y="GeneID",all.y=T)
head(DEG_edgeR_symbol)

# 方法2：gtf文件中得到的id与name关系
# Assembly: GRCh37(hg19) Release: ？
# 使用上课测试得到的count做


# 选择显著差异表达的结果
library(tidyverse)
DEG_edgeR_symbol_Sig <- filter(DEG_edgeR_symbol,regulated!="normal")

# 保存
write.csv(DEG_edgeR_symbol,"DEG_edgeR_triplicate_all.csv", row.names = F)
write.csv(DEG_edgeR_symbol_Sig,"DEG_edgeR_triplicate_Sig.csv", row.names = F)
save(DEG_edgeR_symbol,file = "triplicate_-edgeR_nrDEG.Rdata")


# 加载包
library(pheatmap)
library(tidyverse)

# 加载原始表达矩阵
lname <- load(file = "triplicate_Data.Rdata")
lname

express_cpm1 <- rownames_to_column(as.data.frame(express_cpm) ,var = "ID")

lname <- load(file = "triplicate_-edgeR_nrDEG.Rdata")
lname

# 提取所有差异表达的基因名
edgeR_sigGene <- DEG_edgeR_symbol[DEG_edgeR_symbol$regulated!="normal",]
head(edgeR_sigGene)

data <- merge(edgeR_sigGene,express_cpm1,by.x = "ENSEMBL",by.y = "ID")
data <- na.omit(data)
data <- data[!duplicated(data$SYMBOL),]

# 绘制热图
load(file = "triplicate_Data.Rdata")
load(file = "triplicate_-edgeR_nrDEG.Rdata")
edgeR_sigGene <- DEG_edgeR[DEG_edgeR$regulated!="normal",1]
head(edgeR_sigGene)
dat <- express_cpm[match(edgeR_sigGene,rownames(express_cpm)),]   # 标准化后的表达矩阵，行是基因，列是样本
dat[1:4,1:6]
group <- data.frame(group=group_list)  # 分组信息
rownames(group)=colnames(dat)



ann_colors =list(group = c (A543G = "cyan",WT = "#ff4d40"))
p <- pheatmap(dat,    # 热图数据源，数值型矩阵，numeric
              scale = "row",          # 按行标准化
              show_colnames =T,       # 是否显示列名
              show_rownames = F,      # 是否显示行名
              annotation_colors = ann_colors,
              cluster_cols = T,       # 是否按列聚类
              showannotationname = FALSE,
              annotation_col=group,   # 构建列注释数据框
              main = "Heatmap_of_DEGs")   # 热图标题 
p 

#～～～top100热图p3～～～

#用 https://mp.weixin.qq.com/s/NC5mKwAsYbm1Q-PFC5IEOg 代码
library(pheatmap)

dat <- normal.mtx
FC <- res.df$log2FoldChange
names(FC) <- rownames(res.df)

# 排序差异倍数，提取前100和后100的基因名
DEG_200 <- c(names(head(sort(FC),100)),names(tail(sort(FC),100)))

# 提取基因的归一化
dat <- t(scale(t(dat[DEG_200,])))

group_list <- df_all$group
# 添加注释条
group <- data.frame(group=group_list)
rownames(group) <- colnames(dat)


# 大于2的值赋值为2
dat[dat > 2] <- 2
# 低于-2的值赋值为-2
dat[dat < -2] <- -2
pheatmap(dat, cluster_cols = T,
         show_colnames =F,show_rownames = F,
         annotation_col=group,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))




# ggsave(plot = p3, path = path,'check_top1000_sd.pdf',width = 5,height = 7)
ggsave('qc_top1000_pheatmap.pdf',plot = p3,width = 5,height = 5)




dat <- data[c("WT_1","WT_2","WT_3","A543G_1","A543G_2","A543G_3")]
rownames(dat) <- data$SYMBOL
anno <- data.frame(group=group$sample_title,row.names = group$run_accession)
p <- pheatmap(dat,    # 热图数据源，数值型矩阵，numeric
              scale = "row",          # 按行标准化
              show_colnames =T,       # 是否显示列名
              show_rownames = F,      # 是否显示行名
              cluster_cols = T,       # 是否按列聚类
              annotation_col=group,   # 构建列注释数据框
              main = "edgeR's DEG")   # 热图标题 
p 




lname <- load(file = "triplicate_-edgeR_nrDEG.Rdata")
lname
edgeR_sigGene <- DEG_edgeR_symbol[DEG_edgeR_symbol$regulated!="normal",]
head(edgeR_sigGene)
data <- merge(edgeR_sigGene,express_cpm1,by.x = "ENSEMBL",by.y = "ID")


# 根据需要修改DEG的值
data <- DEG_edgeR_symbol
colnames(data)


# 绘制火山图
colnames(data)


p <- ggplot(data=data, aes(x=logFC, y=-log10(PValue),color=regulated)) + 
  geom_point(alpha=0.5, size=1.5) + 
  theme_set(theme_set(theme_bw(base_size=20))) + theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  xlab("log2FC") + ylab("-log10(Pvalue)") +
  scale_colour_manual(values = c(down='#00AFBB',normal='grey',up='#bb0c00')) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Regulated', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11, hjust = 1)) +
  theme(text = element_text(size = 15)) +
  geom_vline(xintercept=c(-(log2(2)),log2(2)),lty=2,col="black",lwd=0.6) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.6)
p



png(file = "Volcano_Plot.png",width = 900, height = 800, res=150)
plot(p)
dev.off()

#带标签
#创建一个新的列label，并初始化为NA
data$label <- NA

#根据symbol的值，为特定基因添加标签信息
data$label[which(data$ENSEMBL == "ENSDARG00000101244")] <- as.character(expression(italic("gli1")))
data$label[which(data$ENSEMBL == "ENSDARG00000016404")] <- as.character(expression(italic("ptch1")))
data$label[which(data$ENSEMBL == "ENSDARG00000060397")] <- as.character(expression(italic("hhip")))
data$label[which(data$ENSEMBL == "ENSDARG00000062972")] <- as.character(expression(italic("hip1")))
data$label[which(data$ENSEMBL == "ENSDARG00000006837")] <- as.character(expression(italic("mycn")))
data$label[which(data$ENSEMBL == "ENSDARG00000101637")] <- as.character(expression(italic("ccnd1")))
data$label[which(data$ENSEMBL == "ENSDARG00000051748")] <- as.character(expression(italic("ccnd2a")))
data$label[which(data$ENSEMBL == "ENSDARG00000070408")] <- as.character(expression(italic("ccnd2b")))
data$label[which(data$ENSEMBL == "ENSDARG00000052846")] <- as.character(expression(italic("fsta")))

data$label[which(data$ENSEMBL == "ENSDARG00000014246")] <- as.character(expression(italic("jag2a")))
data$label[which(data$ENSEMBL == "ENSDARG00000021389")] <- as.character(expression(italic("jag2b")))
data$label[which(data$ENSEMBL == "ENSDARG00000057992")] <- as.character(expression(italic("fstb")))
data$label[which(data$ENSEMBL == "ENSDARG00000029546")] <- as.character(expression(italic("grem1a")))
data$label[which(data$ENSEMBL == "ENSDARG00000094704")] <- as.character(expression(italic("bcl2a")))
data$label[which(data$ENSEMBL == "ENSDARG00000089109")] <- as.character(expression(italic("bcl2b")))
data$label[which(data$ENSEMBL == "ENSDARG00000055966")] <- as.character(expression(italic("cflara")))
data$label[which(data$ENSEMBL == "ENSDARG00000015399")] <- as.character(expression(italic("foxf1")))
data$label[which(data$ENSEMBL == "ENSDARG00000008133")] <- as.character(expression(italic("foxl1")))
data$label[which(data$ENSEMBL == "ENSDARG00000002445")] <- as.character(expression(italic("prdm1a")))
data$label[which(data$ENSEMBL == "ENSDARG00000098925")] <- as.character(expression(italic("prdm1b")))
data$label[which(data$ENSEMBL == "ENSDARG00000026577")] <- as.character(expression(italic("cdk2")))
data$label[which(data$ENSEMBL == "ENSDARG00000104618")] <- as.character(expression(italic("grem1b")))

data$label[which(data$ENSEMBL == "ENSDARG00000056801")] <- as.character(expression(italic("sufu")))
data$label[which(data$ENSEMBL == "ENSDARG00000033099")] <- as.character(expression(italic("kif7")))
data$label[which(data$ENSEMBL == "ENSDARG00000087538")] <- as.character(expression(italic("kif3a")))
data$label[which(data$ENSEMBL == "ENSDARG00000017174")] <- as.character(expression(italic("dlx2b")))

data$label[which(data$ENSEMBL == "ENSDARG00000013125")] <- as.character(expression(italic("dlx1a")))
data$label[which(data$ENSEMBL == "ENSDARG00000038967")] <- as.character(expression(italic("cul3a")))
data$label[which(data$ENSEMBL == "ENSDARG00000100519")] <- as.character(expression(italic("spop")))
data$label[which(data$ENSEMBL == "ENSDARG00000002952")] <- as.character(expression(italic("smo")))

# 绘制标签火山图
colnames(data)
library(ggrepel)

p1 <- p+ geom_label_repel(data = data, aes(label = label), size = 4, 
                          box.padding = unit (0.5, "lines"),
                          point.padding = unit (0.8, "lines"),
                          segment.color = "black",
                          parse = TRUE,
                          show.legend = FALSE,
                          max.overlaps = 10000)
p1 

png(file = "Volcano_Plot_with_label-2.png",width = 900, height = 800, res=150)
plot(p1)
dev.off()

#--------------------------------------------#
rm(list = ls())
options(stringsAsFactors = F)


library(clusterProfiler)
library(org.Dr.eg.db)
library(GSEABase)
library(ggplot2)
library(tidyverse)

## Error in download.KEGG.Path(species)
# https://github.com/YuLab-SMU/clusterProfiler/pull/471
getOption("clusterProfiler.download.method")

#R.utils::setOption("clusterProfiler.download.method",'auto')
options(clusterProfiler.download.method = "wininet")
#options(clusterProfiler.download.method = "wget")
getOption("clusterProfiler.download.method")

load(file = "triplicate_-edgeR_nrDEG.Rdata")
ls()
# 提取所有差异表达的基因名
DEG <- DEG_edgeR_symbol[DEG_edgeR_symbol$regulated!="normal",2]
head(DEG)

## ===GO数据库, 输出所有结果，后续可根据pvalue挑选结果
ego_CC <- enrichGO(gene=DEG, OrgDb= 'org.Dr.eg.db', keyType='SYMBOL', ont="CC", pvalueCutoff= 1,qvalueCutoff= 1)
ego_MF <- enrichGO(gene=DEG, OrgDb= 'org.Dr.eg.db', keyType='SYMBOL', ont="MF", pvalueCutoff= 1,qvalueCutoff= 1)
ego_BP <- enrichGO(gene=DEG, OrgDb= 'org.Dr.eg.db', keyType='SYMBOL', ont="BP", pvalueCutoff= 1,qvalueCutoff= 1)

p_BP <- barplot(ego_BP,showCategory = 10, label_format = 100) + ggtitle("Biological process")
p_CC <- barplot(ego_CC,showCategory = 10, label_format = 100) + ggtitle("Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10, label_format = 100) + ggtitle("Molecular function")
plotc <- p_BP/p_CC/p_MF
plotc


ggsave('triplicate_enrichGO.png', plotc, width = 10,height = 16)

ego_BP <- data.frame(ego_BP)
ego_CC <- data.frame(ego_CC)
ego_MF <- data.frame(ego_MF)
write.csv(ego_BP,'triplicate_enrichGO_BP.csv')
write.csv(ego_CC,'triplicate_enrichGO_CC.csv')
write.csv(ego_MF,'triplicate_enrichGO_MF.csv')

## === KEGG
genelist <- bitr(gene=DEG, fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Dr.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'zebrafish', pvalueCutoff = 1, qvalueCutoff = 1)
p1 <- barplot(ekegg, showCategory=10, label_format = 100)
p2 <- dotplot(ekegg, showCategory=10, label_format = 100)
plotc = p1/p2
plotc

ggsave('triplicate_enrichKEGG.png', plot = plotc, width = 8, height = 10)

ekegg <- data.frame(ekegg)
write.csv(ekegg,'triplicate_enrichKEGG.csv')

library(patchwork)
p <- (p_BP/p_CC)|(p_MF/p1)
p


ggsave('triplicate_enrichKEGG_GO-2.png',plot = p,width = 16, height = 9)































