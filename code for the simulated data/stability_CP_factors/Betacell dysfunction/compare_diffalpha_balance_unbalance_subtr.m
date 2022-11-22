% This script shows how to compute the factor match score of factors
% extracted from two simulated data sets

% We use Tensor Toolbox

clear all
clc
close all

%% add other apckages to your path
addpath(genpath('.../tensor_toolbox-v3.1')) %Tensor toolbox is needed;  MATLAB Tensor Toolbox. Copyright 2017, Sandia Corporation, http://www.tensortoolbox.org//

%%
load('Fac_CP4_subtr_balance_alpha04.mat','Fac_CP4_perm_subtr')
Fac_CP4_subtr_balance_alpha04=Fac_CP4_perm_subtr;
load('Fac_CP4_subtr_unbalance_alpha04.mat','Fac_CP4_perm_subtr')
Fac_CP4_subtr_unbalance_alpha04=Fac_CP4_perm_subtr;
load('Fac_CP4_subtr_balance_alpha02.mat','Fac_CP4_perm_subtr')
Fac_CP4_subtr_balance_alpha02=Fac_CP4_perm_subtr;
load('Fac_CP4_subtr_unbalance_alpha02.mat','Fac_CP4_perm_subtr')
Fac_CP4_subtr_unbalance_alpha02=Fac_CP4_perm_subtr;


%% compute factor match score--unbalance vs. balance beta=0.2
Fac_CP4_subtr_balance_alpha02_meta_time=ktensor([1;1;1;1],Fac_CP4_subtr_balance_alpha02.U{2},Fac_CP4_subtr_balance_alpha02.U{3});
Fac_CP4_subtr_unbalance_alpha02_meta_time=ktensor([1;1;1;1],Fac_CP4_subtr_unbalance_alpha02.U{2},Fac_CP4_subtr_unbalance_alpha02.U{3});
score_meta_time_unbalance_beta02 = score(Fac_CP4_subtr_unbalance_alpha02_meta_time,Fac_CP4_subtr_balance_alpha02_meta_time,'lambda_penalty',false)
        

%% compute factor match score--balance beta=0.4 vs. balance beta=0.2
Fac_CP4_subtr_balance_alpha04_meta_time=ktensor([1;1;1;1],Fac_CP4_subtr_balance_alpha04.U{2},Fac_CP4_subtr_balance_alpha04.U{3});
score_meta_time_balance_beta04 = score(Fac_CP4_subtr_balance_alpha04_meta_time,Fac_CP4_subtr_balance_alpha02_meta_time,'lambda_penalty',false)


%% compute factor match score--unbalance beta=0.4 vs. balance beta=0.2
Fac_CP4_subtr_unbalance_alpha04_meta_time=ktensor([1;1;1;1],Fac_CP4_subtr_unbalance_alpha04.U{2},Fac_CP4_subtr_unbalance_alpha04.U{3});
score_meta_time_unbalance_beta04 = score(Fac_CP4_subtr_unbalance_alpha04_meta_time,Fac_CP4_subtr_balance_alpha02_meta_time,'lambda_penalty',false)
