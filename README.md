# Infection GWAS

This document describes how GWAS analysis was performed using an example of a single abstract phenotype.

The analysis was performed for EBV, CMV, HSV1, continuous/binary - 6 analyses in total.


## 0. Create config file

In order to prevent putting paths in scripts, we will store the paths to data in a `config` file which looks like this:

```
#!/bin/bash

export geno_dir=''
export files_dir=''

```

This is stored in the base directory.



## 1. Prepare the files for the analysis (phenotypes and covariates)

Notes:

1. The age was selected as the maximum age, for which the measurements were available.

2. Negative values correspond to the missing data.

3. NAs correspond to the withdrawn consent.

4. Binary phenotype files were obtained by applying the threshold 1 for the standardized optical density (dashed red line).


Column names in the raw phenotype data file:

* aln, alnqlet - IDs (pregnancy identifiers): A - 1st child, B - 2nd child, C - 3rd child
* od - optical density
* std - standardized optical density
* std_z - standard z score


```

# cd into the directory with the phenotypic and covariates data files
cd $files_dir

# Load R
module load languages/R-3.5.1-ATLAS-gcc-6.1
R

```

Run

```
Rscript scripts/prepare_phenotypes.r
```

This takes in the data in `files_dir/File_pheno.Rdata` and `files_dir/Gender.txt` and produces both continuous and binary files needed for GWAS analysis for specific type of infection.

## 2. Run GWAS analysis

```
# cd to the directory with the phenotype for the analysis (EBV, CMV, HSV; for each continuous/binary - 6 in total)
cd inf_type/

scripts/run_gwas.sh

```


## 3. Plot the results


#!/bin/bash
```

Rscript scripts/make_plots.r

```




## 4. Clumping

```
scripts/clumping.sh
```

