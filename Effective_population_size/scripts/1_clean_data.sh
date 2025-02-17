plink --bfile ./data/pink_salmon --autosome-num 26 --not-chr 0 --make-bed --out ./work/pink_salmon.clean
plink --bfile ./work/pink_salmon.clean --autosome-num 26 --family --keep-cluster-names Koppen_ODD --hwe .001 --geno 0.02 --maf 0.05 --make-bed --out ./work/Koppen_ODD
plink --bfile ./work/pink_salmon.clean --autosome-num 26 --family --keep-cluster-names Koppen_EVEN --hwe .001 --geno 0.02 --maf 0.05 --make-bed --out ./work/Koppen_EVEN
plink --bfile ./work/pink_salmon.clean --autosome-num 26 --family --keep-cluster-names Nome_ODD --hwe .001 --geno 0.02 --maf 0.05 --make-bed --out ./work/Nome_ODD
plink --bfile ./work/pink_salmon.clean --autosome-num 26 --family --keep-cluster-names Nome_EVEN --hwe .001 --geno 0.02 --maf 0.05 --make-bed --out ./work/Nome_EVEN
plink --bfile ./work/pink_salmon.clean --autosome-num 26 --family --keep-cluster-names Puget_ODD --hwe .001 --geno 0.02 --maf 0.05 --make-bed --out ./work/Puget_ODD
plink --bfile ./work/pink_salmon.clean --autosome-num 26 --family --keep-cluster-names Puget_EVEN --hwe .001 --geno 0.02 --maf 0.05 --make-bed --out ./work/Puget_EVEN
