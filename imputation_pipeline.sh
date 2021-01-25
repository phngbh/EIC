#!/bin/bash

WDIR=/storage/groups/cbm01/workspace/phong.nguyen_new/ethnicity_in_cancer
DDIR=/storage/groups/cbm01/datasets/phong.nguyen/affy_data
REFDIR=/storage/groups/cbm01/datasets/phong.nguyen/reference
LEGFILE=${REFDIR}/legend/1000GP_Phase3_combined.legend
VCFFILE=${REFDIR}/vcf
M3VCFFILE=${REFDIR}/m3vcf
BCFFILE=${REFDIR}/bcf
MAPFILE=${REFDIR}/map/genetic_map_hg19_withX.txt
SNP=${WDIR}/snp.txt

#Set working directory
echo Working directory is ${WDIR} | sed G
cd $WDIR

#Create directories 
echo Creating working folders...

dirlist=(scratch transformed plinkd qc input phased imputed ancestral_snp)
for dir in  ${dirlist[@]} 
do
	if [[ -d ${WDIR}/${dir} ]] 
	then
		echo "Folder ${dir} exists, making new folder..." | sed G
		rm -r ${dir}
		mkdir ${dir}
	else
		echo Making ${dir} | sed G
		mkdir ${dir} 
	fi
done 

#Check reference files
echo Checking reference files...

reffile=(${LEGFILE} ${MAPFILE})
for file in ${reffile[@]}
do
	if [[ -f ${file} ]] 
  then
		echo Reference file ${file} exist | sed G
	else 
		echo "Reference file ${file} does not exist, please provide reference files"
		echo Exiting | sed G
		exit 1
	fi
done

refdir=(${VCFFILE} ${M3VCFFILE} ${BCFFILE})
for dir in ${refdir[@]}
do
	if [[ -d ${dir} ]] 
  then
		echo Reference directory ${dir} exist | sed G
	else 
		echo "Reference directory ${file} does not exist, please make reference directory"
		echo Exiting | sed G
		exit 1
	fi
done

echo Checking input files...

#Check if input files exist 
filename=($(ls $DDIR))
if [[ -z "$filename" ]] 
then
	echo No genotyping file found
	echo Exiting | sed G
	exit 1
 else
   echo Found input files | sed G
fi 

#Check if input files have same length
for file in  ${filename[@]}
do
	nr=$(cat ${DDIR}/${file} | wc -l)
	echo "$file $nr" >> ./scratch/nrow.txt 
done

maxrow=$(cat ./scratch/nrow.txt | awk '{print $2}' | sort -n | tail -n 1)
awk -v max=${maxrow} '$2!=max{print $1}'  ./scratch/nrow.txt > ./scratch/missing_cell.txt 

if [[ -s ./scratch/missing_cell.txt ]] 
then
	echo "Number of SNPs are not equal in all cells. These cells do not have the same number of SNPs:"
	cat ./scratch/missing_cell.txt
	echo Removing these cells... | sed G
	filename=($(cat ./scratch/nrow.txt | awk -v max=${maxrow} '$2==max{print $1}'))
 else
   echo "Number of SNPs are equal in all cells" | sed G
fi 

echo Finish checking | sed G
echo Starting input transformation... | sed G | sed G


for chr in {1..21}
do
	
###TRANSFORMING INPUT AND MAKING PLINK BINARY FILE######
	
  echo Processing chromosome $chr | sed G
	
  echo Checking missing SNPs in chromosome $chr | sed G
	
	for file in  ${filename[@]}
	do 
		nr=$(awk -F "\"*,\"*" -v chr=$chr '$1==chr{print}' ${DDIR}/${file} | wc -l)
		echo "$file $nr" >> ./scratch/nrow_${chr}.txt
	done 
	
	maxrow=$(cat ./scratch/nrow_${chr}.txt | awk '{print $2}' | sort -n | tail -n 1)
	awk -v max=${maxrow} '$2!=max{print $1}'  ./scratch/nrow_${chr}.txt > ./scratch/missing_cell_${chr}.txt 

	if [[ -s ./scratch/missing_cell_${chr}.txt ]] 
	then
		echo "Number of SNPs in chromosome $chr are not equal in all cells. These cells are problematic:"
		cat ./scratch/missing_cell_${chr}.txt
		echo Removing these cells for chromosome $chr ... | sed G
		filename=($(cat ./scratch/nrow_${chr}.txt | awk -v max=${maxrow} '$2==max{print $1}'))
    else
     echo All cells are fine in chromosome ${chr} | sed G
	fi 
	
  echo "Start transforming chromosome ${chr}" | sed G
	
  echo "Making chr${chr}.ped file..." | sed G
  
	for file in  ${filename[@]} 
	do
		#Filter for chromosome and add whitespace within geneotypes
		awk -F "\"*,\"*" -v chr=$chr '$1==chr{print $12}' ${DDIR}/${file} | sed 's/./& /g' > ./scratch/tmp.txt 
		
		#Mark monomorphic and deleted sites with "0"
		awk '{
				if (NF==0)
					print "0 0 ";
				else if (NF==1)
					print "$0 0 ";
				else 
					print $0
				}' ./scratch/tmp.txt > ./scratch/genotype_${file}_${chr}.txt
		
		#Flip column -> row
		paste -s -d "" ./scratch/genotype_${file}_${chr}.txt > ./scratch/genotype_flipped_${file}_${chr}.txt
		echo "$file" > ./scratch/cell_name_tmp.txt
		
		#Remove unneccesary part in cell name
		perl -pi -e 's/_complexGenotypes.csv//g' ./scratch/cell_name_tmp.txt 
		
		#Make .ped file
		paste -d " " ./scratch/cell_name_tmp.txt ./scratch/genotype_flipped_${file}_${chr}.txt >> ./transformed/chr${chr}.ped
   
    rm ./scratch/genotype_flipped_${file}_${chr}.txt
    rm ./scratch/genotype_${file}_${chr}.txt

	done
 
  echo Finish PED file | sed G

	#Make .map file
	echo "Making chr${chr}.map file" | sed G
	awk -F "\"*,\"*" -v chr=$chr '$1==chr{print $1, $6, $2}' ${DDIR}/${filename[0]} > ./transformed/chr"$chr".map
  
  echo "Finished plink text fileset for chromosome ${chr}" | sed G
  
  echo "Making plink binary fileset for chromosome ${chr}" | sed G
	#Make plink binary files
	plink --silent --file ./transformed/chr"$chr" --no-fid --no-parents --no-sex --no-pheno --make-bed --out ./plinkd/chr"$chr"
	
	#Make plink frequency file 
	plink --silent --freq --bfile ./plinkd/chr"$chr" --out ./plinkd/chr"$chr" 
  
  echo "I'm done transforming chromosome ${chr}" | sed G | sed G
  
done 

echo Finish transforming all chromosomes | sed G | sed G 

echo Start quality control | sed G 

###QUALITY CONTROL##########

for chr in {1..21}
do

	echo Start quality control for chromosome ${chr} | sed G

	#Remove SNPs with missing rate >10% 
	echo Checking genotyping rate, MAF and HWE integrity... | sed G
	
	plink --bfile ./plinkd/chr$chr --silent --missing --out ./qc/chr$chr 
 
	sed 1d ./qc/chr${chr}.lmiss | awk '{if ($5 > 0.1) print $2,$5}' > ./qc/chr${chr}_missing_snps.txt
 
	if [[ -s ./qc/chr${chr}_missing_snps.txt ]] 
	then 
		echo "These SNPs have larger than 10% missing rate (see complete list in ./qc/chr${chr}_missing_snps.txt):"
		cat ./qc/chr${chr}_missing_snps.txt | head 
        echo '...' | sed G
	else 
		echo "All SNPs have higher than 90% typing rate" | sed G
	fi
	
	echo "Removing SNPs that have typing rate < 90% ..." | sed G
	plink --bfile ./plinkd/chr${chr} \
    --silent \
    --geno 0.1 --recode --out ./qc/chr${chr}_freqfiltered


	#SNPs with unmatched positions
 
	echo Checking positions compared to reference legend file... | sed G
	
	cat ${LEGFILE} | awk -v CHR=${chr} '{if ($2==CHR && $6=="Biallelic_SNP") print}' > ./qc/ref_chr${chr}.txt
	awk '{print $3}' ./qc/ref_chr${chr}.txt > ./qc/ref_chr${chr}_coord.txt
	echo "Reference legend written to ./qc/ref_chr${chr}.txt and ./qc/ref_chr${chr}_coord.txt" | sed G

	echo Comparing... | sed G

	Rscript --vanilla check_position.R ./qc/chr${chr}_freqfiltered.map ./qc/ref_chr${chr}_coord.txt
	awk '{print $2}' ./qc/unmatched_position.txt > ./qc/chr${chr}_unmatched_position.txt
	
	if [[ -s ./qc/chr${chr}_unmatched_position.txt ]] 
	then 
		echo "SNPs having unmatched positions (see complete list at ./qc/chr${chr}_unmatched_position.txt):"
		cat ./qc/chr${chr}_unmatched_position.txt | head
    	echo '...' | sed G
		echo Removing these SNPs... | sed G
		plink --file ./qc/chr${chr}_freqfiltered --silent --exclude ./qc/chr${chr}_unmatched_position.txt --recode --out ./qc/chr${chr}_posfiltered
	else 
		echo All SNPs have correct positions | sed G
    	plink --file ./qc/chr${chr}_freqfiltered --silent --recode --out ./qc/chr${chr}_posfiltered
	fi 

	#Check REF/ALT and flip strands
  
  echo Checking REF/ALT integrity...
  
  awk '{print $3,$4,$5}' ./qc/ref_chr${chr}.txt > ./qc/ref_chr${chr}_alleles.txt
  echo "Legend alleles are written to ./qc/ref_chr${chr}_alleles.txt" | sed G
  
  echo Checking unmatching and flipping...
  
  rm -f ./qc/swapped_alleles.txt ./qc/unmatched_alleles.txt
  plink --file ./qc/chr${chr}_posfiltered --silent --make-bed --out ./qc/chr${chr}_posfiltered
  Rscript --vanilla check_refalt.R ./qc/chr${chr}_posfiltered.bim ./qc/ref_chr${chr}_alleles.txt
  
  cp ./qc/unmatched_alleles.txt ./qc/chr${chr}_unmatched_alleles.txt
  
  
  if [[ -s ./qc/chr${chr}_unmatched_alleles.txt ]]
  then
    echo "These SNPs have invalid REF/ALT alleles (see complete list at ./qc/chr${chr}_unmatched_alleles.txt):"
    cat ./qc/chr${chr}_unmatched_alleles.txt | head
    echo '...' | sed G
    echo Removing them... | sed G
    plink --bfile ./qc/chr${chr}_posfiltered --silent --exclude ./qc/chr${chr}_unmatched_alleles.txt --make-bed --out ./qc/chr${chr}_QCed
  else 
	  echo No invalid alleles found | sed G
    plink --bfile ./qc/chr${chr}_posfiltered --silent --make-bed --out ./qc/chr${chr}_QCed
  fi
  
  echo "Finish QC check for chromosome ${chr}. Final QCed files at ./qc/chr${chr}_QCed" | sed G | sed G
 
done

echo Finish QC for all chromosomes | sed G | sed G 
echo Start imputation | sed G 

### IMPUTATION #########

echo Creating BCF reference and target files for phasing... | sed G
for chr in {1..21} 
do
	echo Creating indexed BCF files for chromosome ${chr}...
	bcftools view ${VCFFILE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O b -o ${BCFFILE}/ref_chr${chr}.bcf
    bcftools index ${BCFFILE}/ref_chr${chr}.bcf
    plink --silent --bfile ./qc/chr${chr}_QCed --recode vcf-iid bgz --out ./input/chr${chr}
    bcftools view ./input/chr${chr}.vcf.gz -O b -o ./input/chr${chr}.bcf
    bcftools index ./input/chr${chr}.bcf
done

echo Phasing... | sed G 
for chr in {1..21}
do
	echo Phasing chromosome ${chr}...
	eagle --vcfRef ${BCFFILE}/ref_chr${chr}.bcf \
			--vcfTarget ./input/chr${chr}.bcf  \
			--geneticMapFile ${MAPFILE} \
			--outPrefix ./phased/chr${chr}.phased \
			--allowRefAltSwap \
			--numThreads 5 \
			--vcfOutFormat z 
done

echo Imputing...
for chr in {1..21}
do
	echo Imputing chromosome ${chr}...
	minimac4 --refHaps ${M3VCFFILE}/${chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
			 --haps ./phased/chr${chr}.phased.vcf.gz \
             --prefix ./imputed/chr${chr}.imputed \
             --cpus 5  \
             --noPhoneHome \
             --format GT,DS,GP \
             --allTypedSites --minRatio 0.00001
done

echo Finished imputation | sed G | sed G
echo Start post-processing | sed G 

for chr in {1..21}
do
	
	echo Processing imputed chromosome ${chr}
	
	plink --vcf ./imputed/chr${chr}.imputed.dose.vcf.gz \
		--snps-only \
		--silent \
		--recode --out ./scratch/chr${chr}_filtered
	if [[ ! -s ./scratch/chr${chr}_filtered.ped ]]; then 
        echo No snp passing filtering | sed G
		continue 
	fi
  perl -pi -e 's/"//g' $SNP
	awk -v CHR=${chr} '$2==CHR{print}' $SNP > ./scratch/snp_${chr}.txt
	Rscript --vanilla match_ids.R ./scratch/chr${chr}_filtered.map ./scratch/snp_${chr}.txt ${chr}
	
	awk '{print $5}' ./scratch/snp_${chr}_withIDs.txt > ./scratch/extracted_ids_${chr}.txt 
	plink --file ./scratch/chr${chr}_filtered \
		--silent \
		--extract ./scratch/extracted_ids_${chr}.txt \
		--recode vcf-iid --out ./scratch/chr${chr}_filtered_extracted
	awk '!/#|##/ {print}' ./scratch/chr${chr}_filtered_extracted.vcf > ./scratch/chr${chr}_filtered_extracted_trimmed.vcf
	awk '/#CHROM/{print}' ./scratch/chr${chr}_filtered_extracted.vcf > ./ancestral_snp/snp_imputed_${chr}.txt
	#perl -pi -e 's/0_1_//g' ./ancestral_snp/snp_imputed_${chr}.txt
	while IFS= read -r line; do
		echo $line > ./scratch/tmp.txt
		REF=$(awk '{print $4}' ./scratch/tmp.txt)
		ALT=$(awk '{print $5}' ./scratch/tmp.txt)
		perl -pi -e 's/0\/0/'${REF}''${REF}'/g' ./scratch/tmp.txt
		perl -pi -e 's/0\/1/'${REF}''${ALT}'/g' ./scratch/tmp.txt
		perl -pi -e 's/1\/1/'${ALT}''${ALT}'/g' ./scratch/tmp.txt
		cat ./scratch/tmp.txt >> ./ancestral_snp/snp_imputed_${chr}.txt
	done < ./scratch/chr${chr}_filtered_extracted_trimmed.vcf
	perl -pi -e 's/ /\t/g' ./ancestral_snp/snp_imputed_${chr}.txt
done 
