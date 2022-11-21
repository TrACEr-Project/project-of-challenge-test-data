clear all
clc
close all

%% compute factor match score--unbalance vs. balance beta=0.2
%  load factors
load('Fac_CP4_full_balance_alpha02.mat','Fac_CP4_perm_full')
Fac_CP4_full_balance_alpha02=Fac_CP4_perm_full;
load('Fac_CP4_full_unbalance_alpha02.mat','Fac_CP4_perm_full')
Fac_CP4_full_unbalance_alpha02=Fac_CP4_perm_full;
Fac_CP4_full_balance_alpha02_meta_time=ktensor([1;1;1;1],Fac_CP4_full_balance_alpha02.U{2},Fac_CP4_full_balance_alpha02.U{3});
Fac_CP4_full_unbalance_alpha02_meta_time=ktensor([1;1;1;1],Fac_CP4_full_unbalance_alpha02.U{2},Fac_CP4_full_unbalance_alpha02.U{3});
score_meta_time_unbalance_beta02 = score(Fac_CP4_full_unbalance_alpha02_meta_time,Fac_CP4_full_balance_alpha02_meta_time,'lambda_penalty',false)




%% compute factor match score--balance beta=0.4 vs. balance beta=0.2
%  load factors
load('Fac_CP5_full_balance_alpha04.mat','Fac_CP5_perm_full')
Fac_CP5_full_balance_alpha04=Fac_CP5_perm_full;
Fac_CP5_full_balance_alpha04_meta_time=ktensor([1;1;1;1;1],Fac_CP5_full_balance_alpha04.U{2},Fac_CP5_full_balance_alpha04.U{3});
score_meta_time_balance_beta04 = score(Fac_CP5_full_balance_alpha04_meta_time,Fac_CP4_full_balance_alpha02_meta_time,'lambda_penalty',false)

%% compute factor match score--unbalance beta=0.4 vs. balance beta=0.2
% load factors
load('Fac_CP5_full_unbalance_alpha04.mat','Fac_CP5_perm_full')
Fac_CP5_full_unbalance_alpha04=Fac_CP5_perm_full;
Fac_CP5_full_unbalance_alpha04_meta_time=ktensor([1;1;1;1;1],Fac_CP5_full_unbalance_alpha04.U{2},Fac_CP5_full_unbalance_alpha04.U{3});
score_meta_time_unbalance_beta04 = score(Fac_CP5_full_unbalance_alpha04_meta_time,Fac_CP4_full_balance_alpha02_meta_time,'lambda_penalty',false)

