'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-04-19 22:06:17
LastEditors: zpliu
LastEditTime: 2023-04-19 22:17:34
@param: 
'''
import numpy as np
import pandas as pd
import scipy.stats as ss
import sys 
def rank_INT(series, c=0.5, stochastic=False):
    """ Perform rank-based inverse normal transformation on pandas series.
        If stochastic is True ties are given rank randomly, otherwise ties will
        share the same value. NaN values are ignored.
        Args:
            param1 (pandas.Series):   Series of values to transform
            param2 (Optional[float]): Constand parameter (Bloms constant)
            param3 (Optional[bool]):  Whether to randomise rank of ties
        
        Returns:
            pandas.Series
    """

    # Check input
    assert(isinstance(series, pd.Series))
    assert(isinstance(c, float))
    assert(isinstance(stochastic, bool))

    # Set seed
    np.random.seed(123)

    # Take original series indexes
    orig_idx = series.index

    # Drop NaNs
    series = series.loc[~pd.isnull(series)]

    # Get ranks
    if stochastic == True:
        # Shuffle by index
        series = series.loc[np.random.permutation(series.index)]
        # Get rank, ties are determined by their position in the series (hence
        # why we randomised the series)
        rank = ss.rankdata(series, method="ordinal")
    else:
        # Get rank, ties are averaged
        rank = ss.rankdata(series, method="average")

    # Convert numpy array back to series
    rank = pd.Series(rank, index=series.index)

    # Convert rank to normal distribution
    transformed = rank.apply(rank_to_normal, c=c, n=len(rank))
    
    return transformed[orig_idx]

def rank_to_normal(rank, c, n):
    # Standard quantile function
    x = (rank - c) / (n - 2*c + 1)
    return ss.norm.ppf(x) 


if __name__=="__main__":
    residealsFile=sys.argv[1]
    geneNameFile=sys.argv[2]
    geneDataFrame=sys.argv[3]
    outFile=sys.argv[4]

    resideals_data=pd.read_csv("./peer_factor_5/residuals.csv",header=None,index_col=None,sep=",")
    resideals_data_normal=resideals_data.apply(
        lambda x:rank_INT(x),axis=0
    )
    geneId=pd.read_csv("./peer/filter_genePair.txt",header=None,index_col=None,sep="\t")
    SampleList=pd.read_csv("FPKM/gene_expression.txt",header=0,index_col=0,sep="\t")

    resideals_data_normal.index=SampleList.columns
    resideals_data_normal.columns=geneId[0].values
    resideals_data_normal.index.name='IID'
    resideals_data_normal.to_csv(
        "./phenotype/gene_expression.txt",
        header=True,index=True,sep="\t")