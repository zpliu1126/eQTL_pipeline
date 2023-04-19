'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-04-19 22:20:32
LastEditors: zpliu
LastEditTime: 2023-04-19 22:22:32
@param: 
'''
import pandas as pd 
import sys 
#-------------------------------------------------
#TODO: 从gtf文件中制作基因的bed文件
#! 染色体编号需要替换为数字
#-------------------------------------------------

peer_residuals=sys.argv[1]
geneBed=sys.argv[2]
outFile=sys.argv[3]


peer_residuals=pd.read_csv(peer_residuals,header=0,index_col=0,sep="\t")
geneBed=pd.read_csv(geneBed,index_col=4,sep="\t",header=None)
geneBed[0]=geneBed[0].apply(lambda x:int(x.strip("Chr")))
geneBed['TSS']=geneBed.apply(lambda x: x[1] if x[3]=="+" else x[2],axis=1) 


#-------------------------------------------------------------
#TODO 将peer矫正后的数据处理为tensorQTL输入格式
#-------------------------------------------------------------
geneInfo=[]
expressionInfo=[]
for geneID,expressData in peer_residuals.T.iterrows():
    if geneID not in geneBed .index:
        #* 跳过scaffold上的基因
        pass 
    else:
        chrom,start,end,stand,tss=geneBed.loc[geneID]
        geneInfo.append([int(chrom),tss,tss+1,geneID])
        expressionInfo.append(
        expressData.values
        )

geneInfo=pd.DataFrame(geneInfo)
expressionInfo=pd.DataFrame(expressionInfo)
out=pd.concat([geneInfo,expressionInfo],axis=1)
out.columns=['#chr','start','end','phenotype']+list(peer_residuals.T.columns)
#* 按照坐标进行排序
out=out.sort_values(by=['#chr','start'])
out.to_csv(outFile,header=True,index=False,sep="\t")   