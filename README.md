sfon_pqtl

Analyze the Brook Charr genetic map for sex chromosome, heterochiasmy and physiological QTL

### Requirements:   
Note that these are downloaded from within the R script   
R/qtl `http://www.rqtl.org/download/`   
corrplot `https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html`      
rlecuyer
plyr

### Setup:   
Put raw data files into 02_raw_data, including the following:    
genotypes file: `Sfon_female_map_v4.3.loc`
map file: `Sfon_female_map_v4.3.map`   
phenotypes file: `SFQTL_phenotypes-full.qua`  



### Overview:
1. Load map and genotype data, perform quality control and data transformations   
1.1 Subset data to perform heterochiasmy analysis    
1.2 Perform heterochiasmy analysis   
2. Run permutations and create models for QTL analysis
3. Explore and visualize results    

## 1. Load and QC data   
Full commenting is in the R script that indicates all of the steps taken. 
This will result in the saving of the prepared data as: `02_data/sfon_01_output.RData`, which can then be imported into the second script.

## 2. Analyze heterochiasmy   
But before you can run this, you must determine which markers in your dataset are nnxnp markers, and these will be removed in order to only include markers that are efxeg or hkxhk; otherwise there would be a bias as female-specific markers would be present but male specific markers would not.
To do this, run the following command from the terminal:    
`grep -E 'nnxnp' 02_data/Sfon_female_map_v4.3.loc | awk 'BEGIN{ORS=","}1 { print $1 }' | sed 's/,$//g' > 02_data/sfon_nnxnp_markers.csv`
Then you will have a comma-separated list of all of the female-specific (nnxnp) markers to remove from your data.    

Then you can use script `01.1_data_subset*.R` 
This will produce `02_data/sfon_01_output_subset_only_efxeg_and_hkxhk.RData`  

Next, move to the `parentalXO.R` script. Note that this is an adaptation of the plotGeno() (R/qtl) function in order to obtain the parental crossover locations   

Note that there is a special option within this script that lets you run individual chromosomes one-by-one, this was what was used to generate the additional file with individual chromosomes. All that is needed to do to run this version, is to uncomment out the lines in the section `# ##Special (all chromosomes individually)` and comment out the section above it entitled `# # Standard`. Note that the standard version calculates everything separately for the male and the female haplotype, and does the calculation for metacentrics then acrocentrics. This produces the figure seen in the manuscript.
You can also do the same to comment out or access the special version for individual chromosomes in the plotting. 

## 3. Run QTL analysis
Now go back to the dataset with all markers to run the QTL permutations, in script
`02_permutations*.R`
This will give results with chromosome-wide and genome-wide significance, when required using sex as a covariate.
Will produce the following file: `02_data/sfon_02_output.RData`
I suggest you run this on a computer that allows permutations to be parallelized, which can be set in the script. With all of the different phenotypes, a reasonable amount of computing power can speed things up (I ran on MacPro with 18 cores)   

After this is completed, open `03_result_interp*.R` to analyze and visualize your results
Including:   
Results of single QTL, no covariate   
Results of single QTL, consider sex as a covariate   
Results of an interaction effect between sex and QTL to find sex-specific QTL   
Extract results of chromosome-wide significant QTL analysis, which will write to a file `02_data/chr_wide_p0.01_w_10000perms.csv`
Note: to deal with this, I suggest:   
`grep -A1 'lod' 02_data/chr_wide_p0.01_w_10000perms.csv | grep -v '--' - | less`
Results of sex as a binary trait (sex chromosome)   
Effect plots (section still needs some work, but can extract the effect sizes of the different QTL on phenotypes per marker type)   
Determine Percent Variance Explained for each trait, including both chromosome and genome-wide significant QTL    
Effect size estimates per marker type
Some experimental tests to investigate genotype by sex interactions for specific traits   
Identify the actual genotypes   
Plot QTL directly onto genetic map (to produce figure used in paper)   

