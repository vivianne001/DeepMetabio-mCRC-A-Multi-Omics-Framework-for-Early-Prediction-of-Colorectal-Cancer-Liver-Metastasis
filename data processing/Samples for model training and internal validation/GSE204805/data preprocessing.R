
rm(list = ls())
#########################################################################
## GEO数据的下载与存放
rm(list=ls())
library(BiocManager)
library(GEOquery)

# 简单下载GES数据的代码
options(stringsAsFactors = F)
library(limma)
library(GEOquery) # 加载GEOqueryR包
gse="GSE39582" # 选择下载GSE42872
eSet=getGEO(gse,
            destdir='.',
            getGPL=F) # 下载GSE42872，eSet是一个list，包括GSE42872的表达矩阵

# 若以上方法不行，则需在官网下载series_matrix.txt文件
eSet=getGEO(filename = "GSE204805_series_matrix.txt.gz",getGPL = F) # 打开本地下载的GEO文件

# 读入tsv.gz格式文件
# 加载必要的库
library(readr)
# 读取文件
data = read_tsv("GSE204805_merged_hs_mm.tsv.gz")
data[1:10,1:10]
class(data)

data = as.data.frame(data)
colnames(data)

# 删除"H_"前缀
data$Geneid = sub("H_","",data$Geneid)
# 删除"."和"."后面的内容
data$Geneid = sapply(strsplit(data$Geneid, "\\."), "[", 1)

data = data[!duplicated(data$Geneid),]

exp = data

# 读入非matrix.txt.gz格式文件
library(data.table)
# 使用 read.table 读取 .txt.gz 文件
data <- read.table("GSE131418_Consortium_prim_met_GE_probe_level.txt.gz", header = TRUE, sep = "\t")
class(data)
data[1:10,1:10]

data2 = read.table("GSE131418_MCC_prim_met_GE_probe_level.txt.gz",header = T,sep = "\t")
exp = cbind(data,data2)

##
# 获取表达矩阵 （这个数据集的txt.gz没有exp表达数据）
exp=exprs(eSet) # 提取GSE42872的表达矩阵
LL=head(exp) # 查看GSE42872表达矩阵的排列
str(exp) # 查看GSE42872表达矩阵的结构
exp=as.matrix(exp)
boxplot(exp) # 查看各sample是否中位线处于同一水平，若不在，则不在的样本和其他样本存在异常
range(exp)
exp=normalizeBetweenArrays(exp) 
range(exp)

# 判断是否需要Log标准化
qx <- as.numeric(quantile(exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { exp[which(exp <= 0)] <- NaN
exp <- log2(exp+1) }

range(exp) # 检查exp的取值范围

exp=as.data.frame(exp)

rownames(exp) = exp$Geneid
exp = exp[,-1]

# 发现exp有不正常的超额值
# which(exp>500,arr.ind = T) # 查询超额值及其位置
# exp=exp[,-c(42,107)] # 删除具有超额值的病人样本
# boxplot(exp)
#exp[is.na(exp)] = 0 # 将exp中的NA值全部改为0

# range(exp) # range()查看exp的最大/最小值，用以判断exp是否log化

#获取临床矩阵
ph=pData(eSet)#获得临床信息
ph=as.data.frame(ph)
# 分割字符串并选择第一个单词
ph$title = sapply(strsplit(ph$title, ","), "[", 1)
ph$title[1:10]

rownames(ph) = ph$title
ph = ph[,-1]
# -----------------------------------------------------------------------------

## 二、GEO数据分组信息和注释获取
rm(list=ls())
load(file="GSE41258_RAW DATA.Rdata")
exp=as.data.frame(exp)
# 2.注释获得-对应探针ID和基因名称
# 在GPL官网上直接复制-txt-excel
library(readxl)
ID=read_excel("GSE159216 ID.xlsx",sheet = 2,col_names = T) # 
ID=as.data.frame(ID)
ID = ID[-c(1:14),]
colnames(ID) = ID[1,]
ID = ID[-1,]

ID = ID[,c(1,8)]
ID = as.data.frame(ID)
ID[1:3,1:2]

##
# 分割 `gene_assignment` 列中的字符串
split_strings <- strsplit(ID$gene_assignment, " // ")
# 选择每个列表的第二个元素
second_words <- sapply(split_strings, function(x) x[2])
# 将 `gene_assignment` 列替换为第二个单词
ID$gene_assignment <- second_words
# 找出 `gene_assignment` 列中以 "LOC" 开头或等于 "---" 的行
bad_rows <- grepl("^LOC", ID$gene_assignment) | ID$gene_assignment == "---"
# 使用这些行来过滤你的数据框
ID <- ID[!bad_rows, ]

# 使用 na.omit 去除 'GeneSymbol' 列中包含 NA 值的行
ID <- ID[!is.na(ID$gene_assignment), ]

# 3.探针ID-Gene Symbol转化
exp = as.data.frame(exp)
exp$ID=rownames(exp)
exp=merge(exp,ID,by="ID")
exp = as.data.frame(exp)
exp[1:3,280:285]
## 使用 sub 函数来替换 Gene Symbol 列中的 " /// .*" 模式的代码
exp$'Gene Symbol' <- sub(" /// .*", "", exp$'Gene Symbol')

# 找出 `Gene Symbol` 列中的重复值
duplicates <- duplicated(exp$gene_assignment)
# 使用这些重复值来过滤你的数据框
exp <- exp[!duplicates, ]
rownames(exp) = exp$gene_assignment

exp = exp[,-285]
exp = exp[,-1]

### 过滤接受治疗和非人源性样本
ph = ph[ph$`treatment:ch1` %in% "NA",]
ph = ph[ph$characteristics_ch1.1 %in% c("cell line: H"),]

# 统一exp和ph的样本名称
k=intersect(colnames(exp),rownames(ph))
ph=ph[k,]
exp=exp[,k]
identical(colnames(exp),rownames(ph))
print(setequal(colnames(exp), ph$title))

## 4.分组

###
## 2. 设定模型
# 设置任务
invalid_colnames <- rownames(exp)[!make.names(rownames(exp), unique=TRUE) %in% rownames(exp)]
print(invalid_colnames)

# 生成符合R变量命名规则的列名
new_rownames <- make.names(rownames(exp), unique=TRUE)
# 将mpe的列名替换为新的列名
rownames(exp) <- new_rownames
# 标准化exp的名称
rownames(exp) <- gsub("^(.*?)\\..*$", "\\1", rownames(exp))

##
exp$gene_name = rownames(exp)
# 标准化mpe的名称
exp$gene_name <- gsub("^(.*?)\\..*$", "\\1", exp$gene_name)
exp$gene_name[1:100]

exp = exp[!duplicated(exp$gene_name),]
rownames(exp) = exp$gene_name
exp = exp[,-964]

save(exp,ph,file = "GSE204805_Raw_Data.Rdata")
##
grouplist=ifelse(ph$`cell type:ch1`=="LM","metastasis","primary")

# 设置参考水平，确定对照组在前，处理组在后的顺序
grouplist=factor(grouplist,levels = c("metastasis","primary"))
table(grouplist)
##
################################

dat = t(exp)
dat = as.data.frame(dat)
dat$status = grouplist

## 保存
save(grouplist,file = "group_metastasis_primary.Rdata")
save(dat,file = "验证集数据框.Rdata")

