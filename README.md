<!--
 * @Descripttion: 
 * @version: 
 * @Author: zpliu
 * @Date: 2023-04-19 21:56:07
 * @LastEditors: zpliu
 * @LastEditTime: 2023-04-19 22:53:45
 * @@param: 
-->


### 1.原始FPKM文件处理

+ 单个基因表达量大于0.1的accessions 超过5%
+ 单个基因第95%表达水平和5%表达水平相差两倍

> `./FPKM/gene_expression.txt` 为初始的基因表达矩阵，`./peer/filter_gene_expression.txt` 按照表达量过滤后的表达矩阵，
> `./peer/filter_genePair.txt` 过滤后剩下的基因ID

```bash
python filter_expression.py  ./FPKM/gene_expression.txt ./peer/filter_gene_expression.txt ./peer/filter_genePair.txt
```

### 2.peer计算协变量

> 例如按照5个协变量因子进行分析

```bash
module load peer/1.0 
python2.7 peer_interface.py -f peer/filter_gene_expression.txt -n 5 -o peer_factor_5
```

#### 2.1 将peer生成的文件进行标准化

> 注意这里需要使用python3 版本，在module load peer时会自动加载python2，因此需要module unload Python/2.7.15
> 输出文件为 `./phenotype/gene_expression.txt` 即标准化后的表型文件

```bash
mkdir phenotype

python  peer_RINT.py  ./peer_factor_5/residuals.csv ./peer/filter_genePair.txt ./FPKM/gene_expression.txt ./phenotype/gene_expression.txt
```

#### 2.2 根据基因的注释信息向表型文件中添加注释信息

> 由于cis-eQTL分析是基因上下游1Mb的变异进行eQTL分析，因此需要在表型文件中指定基因的TSS信息
> 默认跳过scaffold上的基因，不分析; `./phenotype/gene_expression_peer_5.bed` 文件为过滤后的基因，且包含TSS坐标的信息


```bash
python gene_TSS.py ./phenotype/gene_expression.txt  reference_gene.bed ./phenotype/gene_expression_peer_5.bed 

#* 对文件进行压缩并建索引
module load HTSlib/1.9
bgzip  gene_expression_peer_5.bed 
tabix  -p bed  gene_expression_peer_5.bed.gz
```

### 3.进行cis-eQTL分析

cis eQTL分析中的两种模式：
+ `permutation`找出每个基因最显著的lead SNP
+ `nominal` 对于每个基因，输出其cis区域所有SNP的beta与se值（包括lead SNPs）

>这里脚本使用GPU进行跑会快很多，没有GPU权限可以跟管理员申请 ;
> 由于表型文件是矫正peer协变量后剩余的值，因此跑eQTL分析时，只用PAC的协变量即可

+ `All_eGene_permutate.txt` 表型文件中每个基因lead SNP的显著性信息
+ `Garb_01G000010.cis_qtl_pairs.1.parquet` 该基因的nominal结果

```bash
plinkFile='/cotton/Liuzhenping/diaploidCotton_sQTL/Genotype/filter2_Q600_SNPs_joint_216_number_chr'
phenotype='./phenotype/gene_expression_peer_5.bed.gz'
covariant='./PCA_qcovar.txt'
#* permutate 结果
outFile='All_eGene_permutate.txt'
python cis_mapping.py  ${plinkFile} ${phenotype} ${covariant} ${outFile} p

#* nominal 结果将针对单个基因进行分析，并输出压缩的parquet文件
perfix='Garb_01G000010'
python cis_mapping.py  ${plinkFile} ${phenotype} ${covariant} ${perfix} n 
#* nominal结果 



```

