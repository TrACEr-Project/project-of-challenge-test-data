# Matlab code
This repository contains the descriptions of the datasets and the codes we used in the paper 'Analyzing postprandial metabolomics data using multiway models: A simulation study'.

All the implementations are tested on MacOS version 10.15.3.

## Packages needed for implementation
*  Brett W. Bader, Tamara G. Kolda and others. MATLAB Tensor Toolbox, Version 3.1. Available online at https://www.tensortoolbox.org, 2020
* S. Becker, “L-BFGS-B C code with Matlab wrapper,” 2019. Available online at https://github.com/stephenbeckr/L-BFGS-B-C; see also: R. H. Byrd, P. Lu and J. Nocedal. A Limited Memory Algorithm for Bound Constrained Optimization, (1995), SIAM Journal on Scientific and Statistical Computing , 16, 5, pp. 1190-1208
* C.A. Andersson and R. Bro, 'The N-way Toolbox for MATLAB', Chemometrics & Intelligent Laboratory Systems. 52 (1):1-4, 2000. Available online at http://www.models.life.ku.dk/nwaytoolbox
* Eigenvector Research, DataSet Object, available online at https://eigenvector.com/software/dataset-object/

All the packages should be installed as a subfolder under Matlab path


## Descriptions of the datasets under the folder 'simulated datasets'
*  datasets under the subfolder 'Insulin resistance in Skeletal muscle'
    1. The file named 'Simu_6meta_8time_alpha02_IRM_balance.mat' stores the dataset generated with Insulin resistance in Skeletal muscle as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2
    2. The file named 'Simu_6meta_8time_alpha02_IRM_unbalance.mat' stores the dataset generated with Insulin resistance in Skeletal muscle as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2
    3. The file named 'Simu_6meta_8time_alpha04_IRM_balance.mat' stores the dataset generated with Insulin resistance in Skeletal muscle as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4
    4. The file named 'Simu_6meta_8time_alpha04_IRM_unbalance.mat' stores the dataset generated with Insulin resistance in Skeletal muscle as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4
*  datasets under the subfolder 'Betacell dysfunction' 
    1.   The file named 'Simu_6meta_8time_alpha02_betacell_balance.mat' stores the dataset generated with Beta-cell dysfunction as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2
    2. The file named 'Simu_6meta_8time_alpha02_betacell_unbalance.mat' stores the dataset generated with Beta-cell dysfunction as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2
    3. The file named 'Simu_6meta_8time_alpha04_betacell_balance.mat' stores the dataset generated with Beta-cell dysfunction as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4
    4. The file named 'Simu_6meta_8time_alpha04_betacell_unbalance.mat' stores the dataset generated with Beta-cell dysfunction as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4

## Descriptions of the codes under the folder 'code for the simulated data'

*  The file named 'CP_fulldata.m' is an example code for modeling the full-dynamic (noisy/noiseless) data with the CP model using the tensor toolbox with multiple initialisations

*  The file named 'CP_subtractT0.m' is an example code for modeling the T0-corrected data with CP model using the tensor toolbox with multiple initialisations

*  The file named 'CP_R4_full_split10_check.m' is an example code for checking the replicability of the CP model to the full-dynamic data

*  The file named 'PCA_T0.m' is an example code for  modeling the fasting-state (T0) data with PCA model using the svd function from Matlab

*  codes under the subfolder 'functions'
   1.  The file named 'unique_test_CP.m' is the code for numerically checking the uniqueness of the CP factorization
   2.  The file named 'removesubject.m' is for removing a subset of subjects from the considered dataset
   3.  The file named 'TC.m' is for computing the Tucker congruence

*  codes and data under the subfolder 'stability_CP_factors'
    *  codes and data under the subsubfolder 'Insulin resistance in Skeletal muscle' 
       * The file named 'compare_diffalpha_balance_unbalance_full.m' is for  comparing (computing the factor match scores) the factors extracted from the full-dynamic data from different datasets (low vs. high within-group variation and balanced vs. unbalanced samples) with the between-group variation as the Insulin resistance in Skeletal muscle
       * The file named 'compare_diffalpha_balance_unbalance_subtr.m' is for comparing (computing the factor match score) the factors extracted from the T0-corrected data from different datasets (low vs. high within-group variation and balanced vs. unbalanced samples) with the between-group variation as the Insulin resistance in Skeletal muscle
       * The file named 'Fac_CP4_full_balance_alpha02.mat' stores the factors extracted by the 4-component CP model from the full-dynamic data generated with Insulin resistance in Skeletal muscle as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2
       * The file named 'Fac_CP4_full_balance_alpha04.mat' stores the factors extracted by the 4-component CP model from the full-dynamic data generated with Insulin resistance in Skeletal muscle as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4
       * The file named 'Fac_CP4_full_unbalance_alpha02.mat' stores the factors extracted by the 4-component CP model from the full-dynamic data generated with Insulin resistance in Skeletal muscle as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2
       * The file named 'Fac_CP4_full_unbalance_alpha04.mat' stores the factors extracted by the 4-component CP model from the full-dynamic data generated with Insulin resistance in Skeletal muscle as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4
       * The file named 'Fac_CP4_subtr_balance_alpha02.mat' stores the factors extracted by the 4-component CP model from the T0-corrected data generated with Insulin resistance in Skeletal muscle as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2
       * The file named 'Fac_CP4_subtr_balance_alpha04.mat' stores the factors extracted by the 4-component CP model from the T0-corrected data generated with Insulin resistance in Skeletal muscle as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4
       * The file named 'Fac_CP4_subtr_unbalance_alpha02.mat' stores the factors extracted by the 4-component CP model from the T0-corrected data generated with Insulin resistance in Skeletal muscle as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2
       * The file named 'Fac_CP4_subtr_unbalance_alpha04.mat' stores the factors extracted by the 4-component CP model from the T0-corrected data generated with Insulin resistance in Skeletal muscle as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4
     *  codes and data under the subsubfolder 'Betacell dysfunction'
        * The file named 'compare_diffalpha_balance_unbalance_full.m' is for comparing (computing the factor match scores) the factors extracted from the full-dynamic data from different datasets (low vs. high within-group variation and balanced vs. unbalanced samples) with the between-group variation as the Beta-cell dysfunction
        * The file named 'compare_diffalpha_balance_unbalance_subtr.m' is for comparing (computing the factor match score) the factors extracted from the T0-corrected data from different datasets (low vs. high within-group variation and balanced vs. unbalanced samples) with the between-group variation as the Beta-cell dysfunction
        * The file named 'Fac_CP4_full_balance_alpha02.mat' stores the factors extracted by the 4-component CP model from the full-dynamic data generated with Beta-cell dysfunction as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2
        * The file named 'Fac_CP5_full_balance_alpha04.mat' stores the factors extracted by the 5-component CP model from the full-dynamic data generated with Beta-cell dysfunction as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4
        * The file named 'Fac_CP4_full_unbalance_alpha02.mat' stores the factors extracted by the 4-component CP model from the full-dynamic data generated with Beta-cell dysfunction as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2
        * The file named 'Fac_CP5_full_unbalance_alpha04.mat' stores the factors extracted by the 5-component CP model from the full-dynamic data generated with Beta-cell dysfunction as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4
        * The file named 'Fac_CP4_subtr_balance_alpha02.mat' stores the factors extracted by the 4-component CP model from the T0-corrected data generated with Beta-cell dysfunction as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2
        * The file named 'Fac_CP4_subtr_balance_alpha04.mat' stores the factors extracted by the 4-component CP model from the T0-corrected data generated with Beta-cell dysfunction as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4
        * The file named 'Fac_CP4_subtr_unbalance_alpha02.mat' stores the factors extracted by the 4-component CP model from the T0-corrected data generated with Beta-cell dysfunction as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2
        * The file named 'Fac_CP4_subtr_unbalance_alpha04.mat' stores the factors extracted by the 4-component CP model from the T0-corrected data generated with Beta-cell dysfunction as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4

## Descriptions of the codes under the folder 'code for the real data'

*  The file named 'CP_R3_sixselected_metabolites.m' is the code for modeling the T0-corrected real data (with six selected metabolites) with CP model using the tensor toolbox with multiple initialisations
*  The file named 'CP_R3_sixselected_metabolites_split10_check.m' is an example code for checking the replicability of the CP model to the T0-corrected real data


* codes under the subfolder 'functions'
   *  The file named 'unique_test_CP.m' is an example code for numerically checking the uniqueness of the CP factorization
   * The file named 'removesubject.m' is for removing a subset of subjects from the considered dataset
   * The file named 'TC.m' is for computing the Tucker congruence
   * The file named 'removeisnan.m' is used to remove the subjects and features with 70% missing values
   * codes under the subsubfolder 'spca_wopt_functions'
    
**Original repository:** https://github.com/Lu-source/project-of-challenge-test-data
