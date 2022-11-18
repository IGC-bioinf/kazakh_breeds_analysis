plink --bfile ${1} --missing --chr-set 35 --double-id
Rscript --no-save hist_miss.R 
plink --bfile ${1} --geno 0.1 --make-bed --out ${1}_2 -chr-set 35 --double-id
plink --bfile ${1}_2 --mind 0.1 --make-bed --out ${1}_3 -chr-set 35 --double-id
plink --bfile ${1}_3 --geno 0.02 --make-bed --out ${1}_4 -chr-set 35 --double-id
plink --bfile ${1}_4 --mind 0.02 --make-bed --out ${1}_5 -chr-set 35 --double-id
Rscript --no-save gender_check.R
plink --bfile ${1}_5 --freq --out MAF_check -chr-set 35 --double-id
plink --bfile ${1}_5 --maf 0.05 --make-bed --out ${1}_8 -chr-set 35 --double-id
plink --bfile ${1}_8 --hwe 1e-6 --hwe-all --make-bed --out ${1}_9 -chr-set 35 --double-id

#final_step
plink --bfile ${1}_9 --assoc --out assoc_results -chr-set 35 --double-id


Rscript --no-save Manhattan_plot.R
Rscript --no-save QQ_plot.R






