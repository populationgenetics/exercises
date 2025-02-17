# ./scripts/2_do_PCA.sh 
All commands to be run in the terminal.

For each of the two basic data files (before and after removing loci not placed on a chromosome), use Plink to perform a principal components analysis.

```bash

plink --bfile ./data/pink_salmon --autosome-num 26 --maf 0.1 --pca 3 --out ./work/pink_data.initial

plink --bfile ./work/pink_salmon.clean --autosome-num 26 --maf 0.1 --pca 3 --out ./work/pink_salmon.clean

```
