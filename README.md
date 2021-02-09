# Work flow of the EIC project data analysis

1. Create working directory, data directory and reference genome directory with structure as below as below:

          |--- working dir
          |--- data dir
          |--- reference dir 
               |--- legend/1000GP_Phase3_combined.legend
               |--- vcf
               |--- m3vcf
               |--- bcf
               |--- map/genetic_map_hg19_withX.txt

2. Download 1000 Genome Project's VCF files for all populations (ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) and put to ref dir/vcf/ . Download legend file (https://www.well.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz) and put to ref dir/legend/ . Download genetic map file (https://alkesgroup.broadinstitute.org/Eagle/downloads/tables) and put to ref dir/map. 

3. Put 100spns.txt (ancestral SNPs) as well as all the scripts in this repository to the working directory

4. Install the following software in advance: R, R packages ggplot2, dplyr, biomaRt and stringr, Perl, Ensembl API (https://m.ensembl.org/info/docs/api/index.html), Eagle2, minimac3, minimac4

5. Run get_1000Gfreg.pl to get genotype frequencies of the ancestral SNPs using the Ensembl API

6. Run get_snp_coordinates.R to get the SNP coordinates

7. Example command for the imputation pipeline using a docker image containing all the required software (see sample Dockerfile above) and running on a HPC:

              bash submit_job.sh <Docker image> imputation_pipeline.sh <working directory> <data directory> <reference directory>

8. Run ancestry_inference.R to get the final ancestry matrix 
    
       
