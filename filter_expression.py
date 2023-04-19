'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-04-19 21:57:09
LastEditors: zpliu
LastEditTime: 2023-04-19 21:58:36
@param: 
'''
import pandas as pd 
#----------------------------------------------------------------------
#? 基因表达量文件跑peer时，不需要header与index
#* 对基因表达量进行进一步过滤
#? 表达量大于0.1 的样本占比超过5% 
#? 存在两倍表达差异
#----------------------------------------------------------------------
geneExpression=pd.read_csv("FPKM/gene_expression.txt",header=0,index_col=0,sep="\t")
geneExpression=geneExpression.T
lowIndex=int(geneExpression.shape[0]*0.05)
highIndex=int(geneExpression.shape[0]*0.95)
fiterGenePair=[]
for gene,ratio in geneExpression.apply(lambda x: len([i for i in x if i>0.1])/geneExpression.shape[0],axis=0).iteritems():
        if ratio>=0.05:
            lowExpression=geneExpression[gene].sort_values()[lowIndex]
            highExpression=geneExpression[gene].sort_values()[highIndex]
            if (highExpression-lowExpression)/(highExpression+lowExpression)>=0.33:
                #* 95%和5%的样本存在表达差异
                fiterGenePair.append(gene) 

fiter_expression=geneExpression[fiterGenePair]
fiter_expression.to_csv("./peer/filter_gene_expression.txt",header=False,index=False,sep=",")
pd.Series(fiterGenePair).to_csv("peer/filter_genePair.txt",header=False,index=False,sep="\t")