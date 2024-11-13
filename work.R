#肿瘤通路富集分析脚本
#环境配置，本人使用R4.4.1版本，首先使用getwd()命令获得工作目录，使用setwd("path")设置工作目录，将文件拷贝至该目录
#包安装好的不需要再装
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")  
BiocManager::install(version = "3.19")
BiocManager::install("EnhancedVolcano")
BiocManager::install("clusterProfiler")  
BiocManager::install("org.Hs.eg.db")  
BiocManager::install("pathview")
BiocManager::install("DESeq2")
BiocManager::install("pheatmap")
library('DESeq2')
library('EnhancedVolcano')
setwd('/Users/zhang/Library/Mobile Documents/com~apple~CloudDocs/工作文件/生物信息学课程作业')
#数据清洗
countData <- read.table(file = 'geneCountMatrix.txt',header = T ,row.names = 1)
colData <- read.csv('samplesinfo.csv')
#数据处理
dds <- DESeqDataSetFromMatrix(countData = countData,  
                              colData = colData,  
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "PTC", "ATC"), alpha = 0.05)  
# 查看结果,本人选择基因KRAS是第9185位  
print(res[9185, ])
# 导出结果到CSV文件  
write.csv(as.data.frame(res), "DESeq2_results.csv")

#火山图绘图
p1 <- EnhancedVolcano(res,  
                lab = rownames(res), # 使用基因名称作为标签  
                x = 'log2FoldChange', # x轴是log2 fold change  
                y = 'pvalue', # y轴是p值（你可以使用-log10(pvalue)来获得更好的可视化效果）  
                pCutoff = 0.05, # 显著性水平的阈值  
                FCcutoff = 1, # log2 fold change的阈值，确定哪些点将被高亮显示  
                pointSize = 2.0, # 点的大小  
                labSize = 3.0, # 标签的大小  
                col = c('grey', 'red', 'blue'), # 点的颜色：不显著、上调、下调  
                colAlpha = 0.8, # 点的透明度  
                legendPosition = 'right', # 图例的位置  
                legendLabSize = 10, # 图例文字的大小  
                legendIconSize = 4.0, # 图例图标的大小  
                drawConnectors = TRUE, # 是否绘制连接标签和点的线  
                widthConnectors = 0.5, # 连接线的宽度  
                colConnectors = 'grey') # 连接线的颜色  
ggsave(filename = "volcanoplot.png", plot = p1, width = 8, height = 6)# 保存图片
library(clusterProfiler)  
library(org.Hs.eg.db)  # 对于人类基因  
library(pathview)
# 获取显著差异表达的基因（假设阈值为log2FoldChange > 2且padj < 0.05）  
sigGenes <- res[which(res$log2FoldChange > 2 & res$padj < 0.05), ]  
geneList <- rownames(sigGenes)  # 获取这些基因的基因名
ego <- enrichGO(gene          = geneList,  
                OrgDb         = org.Hs.eg.db,  
                keyType       = 'SYMBOL',  
                ont           = "BP",  
                pAdjustMethod = "BH",  
                pvalueCutoff  = 0.05,  
                qvalueCutoff  = 0.2)

geneEntrezIDs <- bitr(geneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID  
ekegg <- enrichKEGG(gene         = geneEntrezIDs,  
                    organism     = 'hsa',  
                    pvalueCutoff = 0.05)

#go富集分析绘图
gob_bar <- barplot(ego, showCategory=10, title="GO Biological Process Enrichment")
gob_dot <- dotplot(ego, showCategory=20, title="GO Biological Process Enrichment")
ggsave(filename = "go_enrichment_bar.png", plot = gob_bar, width = 8, height = 6)
ggsave(filename = "go_enrichment_dot.png", plot = gob_dot, width = 12, height = 9)

dds.res.FC <- sigGenes[, "log2FoldChange"]
names(dds.res.FC) <- rownames(sigGenes)
gob_cnet <- cnetplot(ego, foldChange = dds.res.FC)
ggsave(filename = "go_enrichment_cnetplot.png", plot = gob_cnet, width = 12, height = 9)

#kegg富集分析
kegg_bar <- barplot(ekegg, showCategory=10, title="KEGG Pathway Enrichment")
kegg_dot <- dotplot(ekegg, showCategory=20, title="KEGG Pathway Enrichment")
ggsave(filename = "kegg_enrichment_bar.png", plot = kegg_bar, width = 8, height = 6)
ggsave(filename = "kegg_enrichment_dot.png", plot = kegg_dot, width = 12, height = 9)

dds.res.FC <- sigGenes[, "log2FoldChange"]
names(dds.res.FC) <- rownames(sigGenes)
kegg_cnet <- cnetplot(ekegg, foldChange = dds.res.FC)
ggsave(filename = "kegg_enrichment_cnetplot.png", plot = kegg_cnet, width = 12, height = 9)

#作者：张奕博保留所有权利
