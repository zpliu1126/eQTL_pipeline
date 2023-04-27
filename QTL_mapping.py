'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-12-05 21:29:41
LastEditors: zpliu
LastEditTime: 2023-04-27 22:57:56
@param: 
'''
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans
import os
import torch
import sys

os.environ['CUDA_VISIBLE_DEVICES'] = "0"
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print('using {}'.format(device))

plink_prefix_path = sys.argv[1]
expression_bed = sys.argv[2]
covariates_file = sys.argv[3]
OutFile = sys.argv[4]

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)

covariates_df = pd.read_csv(
    covariates_file, sep='\t', index_col=0, header=None)  # samples x covariates
# * 协变量文件，最终格式为 samples x covariates
covariates_df.drop([1], axis=1)
covariates_df = covariates_df.drop([1], axis=1)
covariates_df.index.name = 'id'
covariates_df.columns = ['PC1', 'PC2', 'PC3']
covariates_df = covariates_df.loc[phenotype_df.columns]

#* 基因型文件, 对应时期的样本编号
pr = genotypeio.PlinkReader(
    plink_prefix_path, select_samples=phenotype_df.columns)

# load genotypes and variants into data frames
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

mode = sys.argv[5]
if mode == 'p':
    # 进行cis-sQTL分析
    cis_df = cis.map_cis(genotype_df,
                         variant_df,
                         phenotype_df,
                         phenotype_pos_df,
                         covariates_df,
                         nperm=1000,
                         maf_threshold=0.05,
                         window=1000000,
                         seed=2022
                         )
    cis_df.to_csv(OutFile, header=True, index=True, sep="\t")
elif mode == 'n':
    geneList=pd.read_csv(sys.argv[6],header=None,index_col=None,sep="\t")
    try:
        filter_phenotype=phenotype_df.loc[geneList[0].values]
        filter_pos=phenotype_pos_df.loc[geneList[0].values]
    except KeyError:
        print("gene ID not in the phenotype!!")
    cis.map_nominal(genotype_df, variant_df,
                             filter_phenotype,
                             filter_pos,
                             prefix=OutFile,
                             covariates_df=covariates_df,
                             maf_threshold=0.05,
                             window=1000000,
                             output_dir='.',
                             write_top=True, write_stats=True)
elif mode=='t':
    #! 保留所有位点的beta与se值
    geneList=pd.read_csv(sys.argv[6],header=None,index_col=None,sep="\t")
    try:
        filter_phenotype=phenotype_df.loc[geneList[0].values]
    except KeyError:
        print("gene ID not in the phenotype!!")
    trans_df = trans.map_trans(genotype_df, filter_phenotype, covariates_df,
                            return_sparse=True, pval_threshold=1, maf_threshold=0.05,
                            batch_size=20000)
    trans_df.to_csv(OutFile,header=True,index=False)


