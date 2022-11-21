% This file lists the descriptions of the datasets and the codes we used in
% the paper ''Analyzing postprandial metabolomics data using multiway models: A simulation study''


%% Descriptions of the datasets under the folder 'simulated datasets'

%%%%%  codes and data under the subfolder 'Insulin resistance in Skeletal muscle'
%  The file named 'Simu_6meta_8time_alpha02_IRM_balance.mat' stores the dataset generated with Insulin resistance in Skeletal muscle
%  as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2. 

%  The file named 'Simu_6meta_8time_alpha02_IRM_unbalance.mat' stores the dataset generated with Insulin resistance in Skeletal muscle
%  as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2. 

%  The file named 'Simu_6meta_8time_alpha04_IRM_balance.mat' stores the dataset generated with Insulin resistance in Skeletal muscle
%  as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4. 

%  The file named 'Simu_6meta_8time_alpha04_IRM_unbalance.mat' stores the dataset generated with Insulin resistance in Skeletal muscle
%  as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4. 


%%%%%  codes and data under the subfolder 'Betacell dysfunction' 
%  The file named 'Simu_6meta_8time_alpha02_betacell_balance.mat' stores the dataset generated with Beta-cell dysfunction 
%  as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2. 

%  The file named 'Simu_6meta_8time_alpha02_betacell_unbalance.mat' stores the dataset generated with Beta-cell dysfunction
%  as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2. 

%  The file named 'Simu_6meta_8time_alpha04_betacell_balance.mat' stores the dataset generated with Beta-cell dysfunction
%  as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4. 

%  The file named 'Simu_6meta_8time_alpha04_betacell_unbalance.mat' stores the dataset generated with Beta-cell dysfunction
%  as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4. 




%% Descriptions of the codes under the folder 'code for the simulated data'

%  The file named 'CP_unique_fulldata.m' is an example code for modeling the full-dynamic data
%  with CP model using the tensor toolbox with multiple initialisations.
%  To run the code, tensor toolbox as well as
%  the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C is needed.

%  The file named 'CP_unique_subtractT0.m' is an example code for modeling the T0-corrected data
%  with CP model using the tensor toolbox with multiple initialisations.
%  To run the code, tensor toolbox as well as
%  the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C is needed.


%  The file named 'CP_R4_full_split10_check.m' is an example code for
%  checking the replicability of the CP model to the full-dynamic data.
%  To run the code, tensor toolbox as well as
%  the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C is needed.
%  The function removesubject.m is also needed to remove a subset of subjects from the considered dataset.


%  The file named 'PCA_T0.m' is an example code for  modeling the fasting-state (T0) data
%  with PCA model using the svd function from Matlab.


%%%%%  codes and data under the subfolder 'Insulin resistance in Skeletal muscle'
%  The file named 'compare_diffalpha_balance_unbalance_full.m' is for
%  comparing (computing the factor match scores) the factors extracted from
%  the full-dynamic data from different datasets (low vs. high within-group variation and balanced
%  vs. unbalanced samples) with the between-group variation as the Insulin
%  resistance in Skeletal muscle.
%  To run the code, tensor toolbox is needed.

%  The file named 'compare_diffalpha_balance_unbalance_subtr.m' is for
%  comparing (computing the factor match score) the factors extracted from
%  the T0-corrected data from different datasets (low vs. high within-group variation and balanced
%  vs. unbalanced samples) with the between-group variation as the Insulin
%  resistance in Skeletal muscle.
%  To run the code, tensor toolbox is needed.

% The file named 'Fac_CP4_full_balance_alpha02.mat' stores the factors
% extracted by the 4-component CP model from the full-dynamic data generated with Insulin resistance in Skeletal muscle
% as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2. 

% The file named 'Fac_CP4_full_balance_alpha04.mat' stores the factors
% extracted by the 4-component CP model from the full-dynamic data generated with Insulin resistance in Skeletal muscle
% as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4. 

% The file named 'Fac_CP4_full_unbalance_alpha02.mat' stores the factors
% extracted by the 4-component CP model from the full-dynamic data generated with Insulin resistance in Skeletal muscle
% as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2. 

% The file named 'Fac_CP4_full_unbalance_alpha04.mat' stores the factors
% extracted by the 4-component CP model from the full-dynamic data generated with Insulin resistance in Skeletal muscle
% as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4. 

% The file named 'Fac_CP4_subtr_balance_alpha02.mat' stores the factors
% extracted by the 4-component CP model from the T0-corrected data generated with Insulin resistance in Skeletal muscle
% as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2. 

% The file named 'Fac_CP4_subtr_balance_alpha04.mat' stores the factors
% extracted by the 4-component CP model from the T0-corrected data generated with Insulin resistance in Skeletal muscle
% as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4. 

% The file named 'Fac_CP4_subtr_unbalance_alpha02.mat' stores the factors
% extracted by the 4-component CP model from the T0-corrected data generated with Insulin resistance in Skeletal muscle
% as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2. 

% The file named 'Fac_CP4_subtr_unbalance_alpha04.mat' stores the factors
% extracted by the 4-component CP model from the T0-corrected data generated with Insulin resistance in Skeletal muscle
% as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4. 




%%%%%  codes and data under the subfolder 'Betacell dysfunction'
%  The file named 'compare_diffalpha_balance_unbalance_full.m' is for
%  comparing (computing the factor match scores) the factors extracted from
%  the full-dynamic data from different datasets (low vs. high within-group variation and balanced
%  vs. unbalanced samples) with the between-group variation as the Beta-cell dysfunction.
%  To run the code, tensor toolbox is needed.

%  The file named 'compare_diffalpha_balance_unbalance_subtr.m' is for
%  comparing (computing the factor match score) the factors extracted from
%  the T0-corrected data from different datasets (low vs. high within-group variation and balanced
%  vs. unbalanced samples) with the between-group variation as the Beta-cell dysfunction.
%  To run the code, tensor toolbox is needed.

% The file named 'Fac_CP4_full_balance_alpha02.mat' stores the factors
% extracted by the 4-component CP model from the full-dynamic data generated with Beta-cell dysfunction
% as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2. 

% The file named 'Fac_CP5_full_balance_alpha04.mat' stores the factors
% extracted by the 5-component CP model from the full-dynamic data generated with Beta-cell dysfunction
% as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4. 

% The file named 'Fac_CP4_full_unbalance_alpha02.mat' stores the factors
% extracted by the 4-component CP model from the full-dynamic data generated with Beta-cell dysfunction
% as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2. 

% The file named 'Fac_CP5_full_unbalance_alpha04.mat' stores the factors
% extracted by the 5-component CP model from the full-dynamic data generated with Beta-cell dysfunction
% as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4. 

% The file named 'Fac_CP4_subtr_balance_alpha02.mat' stores the factors
% extracted by the 4-component CP model from the T0-corrected data generated with Beta-cell dysfunction
% as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.2. 

% The file named 'Fac_CP4_subtr_balance_alpha04.mat' stores the factors
% extracted by the 4-component CP model from the T0-corrected data generated with Beta-cell dysfunction
% as the between-group variation (50 control and 50 diseased subjects) and the within-group variation at level alpha=0.4. 

% The file named 'Fac_CP4_subtr_unbalance_alpha02.mat' stores the factors
% extracted by the 4-component CP model from the T0-corrected data generated with Beta-cell dysfunction
% as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.2. 

% The file named 'Fac_CP4_subtr_unbalance_alpha04.mat' stores the factors
% extracted by the 4-component CP model from the T0-corrected data generated with Beta-cell dysfunction
% as the between-group variation (70 control and 30 diseased subjects) and the within-group variation at level alpha=0.4. 



%% Descriptions of the codes under the folder 'code for the real and simulated data'

%  The file named 'unique_test_CP.m' is an example code for numerically checking the uniqueness of the CP factorization.
%  To run the code, tensor toolbox is needed.

%  The file named 'removesubject.m' is for removing a subset of subjects from the considered dataset.

%  The file named 'TC.m' is for computing the Tucker congruency




%% Descriptions of the codes under the folder 'code for the real data'

%  The file named 'CP_T0corrected_LargeMetaboliteset_Males.m' is the
%  code for  modeling the T0-corrected real data (only male subjects and 162 measurements)
%  with CP model using the tensor toolbox with multiple initialisations.
%  To run the code, tensor toolbox as well as
%  the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C is needed.
%  To run the code, the function removesubject.m is used to remove the outlier subjects, 
%  the function unique_test_CP.m is used to numerically check the uniqueness of the CP factorization, 
%  the function removeisnan.m is used to remove the subjects with lots of missing values, and the function TC.m is used to
%  compute the Tucker congruence. 


%  The file named 'CP_R2_T0subtr_split10_check.m' is an example code for
%  checking the replicability of the CP model to the T0-corrected real data.
%  To run the code, tensor toolbox as well as
%  the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C is needed.
%  The function removesubject.m is also needed to remove a subset of
%  subjects from the considered dataset.

%  The file named 'spca_wopt_T0_LargeMetaboliteset_Males_Missvalues_approxiamtion' is for
%  computing the approximations of the missing values in the fasting-state
%  real data using the weighted optimization (see the work ' Scalable tensor factorizations for incomplete data' by E. Acar, 
%  D.M. Dunlavy, T.G. Kolda, and M. Mørup).
%  To run the code, the Poblano Toolbox from https://github.com/sandialabs/poblano_toolbox is needed. 
%  To run the code, the functions in the folder 'spca_wopt functions' are
%  needed, and the function removesubject.m is used to remove the outlier
%  subjects.


%  The file named 'PCA_T0_LargeMetaboliteset_Males_Missvalues_replaced.m'
%  is for modeling the fasting-state real data with PCA model using the svd function from Matlab, 
%  missing values being replaced by approximations.
%  To run the code, the function removesubject.m is used to remove the outlier subjects.




