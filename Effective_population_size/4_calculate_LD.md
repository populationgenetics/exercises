# ./scripts/4_calculate_LD.sh 
All commands to be run in the terminal.

See Plink documentation [here](https://www.cog-genomics.org/plink2/).

For each of the six poulations, use Plink to calculate LD between each pair of loci.   The "--r2 square" flag to Plink will produce a LxL matrix of pairwise r<sup>2</sup>, with L as the number of loci.  



```bash
plink --bfile ./work/Koppen_ODD --r2 square --out ./work/Koppen_ODD

plink --bfile ./work/Koppen_EVEN --r2 square --out ./work/Koppen_EVEN

plink --bfile ./work/Nome_ODD --r2 square --out ./work/Nome_ODD

plink --bfile ./work/Nome_EVEN --r2 square --out ./work/Nome_EVEN

plink --bfile ./work/Puget_ODD --r2 square --out ./work/Puget_ODD

plink --bfile ./work/Puget_EVEN --r2 square --out ./work/Puget_EVEN

```