
clear all
clc
close all
%%
load('Fac_CP4_perm_subtr_balance_beta04.mat','Fac_CP4_perm_subtr')
Fac_CP4_perm_subtr_balance_beta04=Fac_CP4_perm_subtr;
load('Fac_CP4_perm_subtr_unbalance_beta04.mat','Fac_CP4_perm_subtr')
Fac_CP4_perm_subtr_unbalance_beta04=Fac_CP4_perm_subtr;
load('Fac_CP4_perm_subtr_balance_beta02.mat','Fac_CP4_perm_subtr')
Fac_CP4_perm_subtr_balance_beta02=Fac_CP4_perm_subtr;
load('Fac_CP4_perm_subtr_unbalance_beta02.mat','Fac_CP4_perm_subtr')
Fac_CP4_perm_subtr_unbalance_beta02=Fac_CP4_perm_subtr;


%% compute factor match score--unbalance vs. balance beta=0.2
Fac_CP4_perm_subtr_balance_beta02_meta_time=ktensor([1;1;1;1],Fac_CP4_perm_subtr_balance_beta02.U{2},Fac_CP4_perm_subtr_balance_beta02.U{3});
Fac_CP4_perm_subtr_unbalance_beta02_meta_time=ktensor([1;1;1;1],Fac_CP4_perm_subtr_unbalance_beta02.U{2},Fac_CP4_perm_subtr_unbalance_beta02.U{3});
score_meta_time_unbalance_beta02 = score(Fac_CP4_perm_subtr_unbalance_beta02_meta_time,Fac_CP4_perm_subtr_balance_beta02_meta_time,'lambda_penalty',false)
        

%% compute factor match score--balance beta=0.4 vs. balance beta=0.2
Fac_CP4_perm_subtr_balance_beta04_meta_time=ktensor([1;1;1;1],Fac_CP4_perm_subtr_balance_beta04.U{2},Fac_CP4_perm_subtr_balance_beta04.U{3});
score_meta_time_balance_beta04 = score(Fac_CP4_perm_subtr_balance_beta04_meta_time,Fac_CP4_perm_subtr_balance_beta02_meta_time,'lambda_penalty',false)


%% compute factor match score--unbalance beta=0.4 vs. balance beta=0.2
Fac_CP4_perm_subtr_unbalance_beta04_meta_time=ktensor([1;1;1;1],Fac_CP4_perm_subtr_unbalance_beta04.U{2},Fac_CP4_perm_subtr_unbalance_beta04.U{3});
score_meta_time_unbalance_beta04 = score(Fac_CP4_perm_subtr_unbalance_beta04_meta_time,Fac_CP4_perm_subtr_balance_beta02_meta_time,'lambda_penalty',false)

