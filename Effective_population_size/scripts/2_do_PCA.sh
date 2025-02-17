plink --bfile ./data/pink_salmon --autosome-num 26 --maf 0.1 --pca 3 --out ./work/pink_data.initial
plink --bfile ./work/pink_salmon.clean --autosome-num 26 --maf 0.1 --pca 3 --out ./work/pink_salmon.clean
