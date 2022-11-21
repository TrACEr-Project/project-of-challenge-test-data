% This script shows how to use weighted optimization for fitting PCA to incomplete data 
% (the fasting-state real data with only male subjects and 162 measurements), similarly as the way to fit the CP model to incomplete data 
% introduced in the paper ' Scalable tensor factorizations for incomplete data' by E. Acar, 
% D.M. Dunlavy, T.G. Kolda, and M. Mørup.   

% Poblano Toolbox from https://github.com/sandialabs/poblano_toolbox is used to solve the optimization problem. 
% In addition, parts of the scripts may require the dataset object (https://eigenvector.com/software/dataset-object/), publically available.





clear all
clc
close all
%%
folder = fileparts(which('spca_wopt_T0_LargeMetaboliteset_Males_Missvalues_approxiamtion')); 
addpath(genpath(folder));
load('NMR+Gene+SNP+Metav3_July4_2022.mat') %% load data

%% remove outliers
pid_list = str2num(NMR.label{1});
outlier_pid =[142,79,342,250,90,335,312];
NMR_remove=remove_outliers(NMR,outlier_pid);

%% remove the metabolites with lots of missing values
NMR_remove=removeisnan(NMR_remove);

%% select the large set of metabolites: 1
selected_metabolite_index = find(NMR_remove.class{3,1}==1 | NMR_remove.class{3,1}==2);
NMR_remove= NMR_remove(:,:,selected_metabolite_index);

%% choose male(2)/female(1)
selected_gender_index =find(NMR_remove.class{1,1}==2);
NMR_remove= NMR_remove(selected_gender_index,:,:);


%% fasting-state data (raw)
X_slect = NMR_remove.data(:,1,:);

%% preprocessing
s=size(X_slect);
X_center=X_slect-repmat(nanmean(X_slect,1),s(1),1);
for i=1:s(3)
    temp=X_center(:,:,i);
    X_pre(:,:,i)=temp/sqrt(nanmean(temp.^2,'all'));
end
X=X_pre;
X=permute(X, [1 3 2]);
X=tensor(X);

%% fit spca_wopt function
% create missing-value index matrix
W=ones(size(X));
W(find(isnan(X.data)))=0;
W=tensor(W);
X(find(isnan(X.data)))=0;
% set parameters for spca_wopt 
nb_starts =60;
nm_comp=2;
options=ncg('defaults');
options.DisplayIters=1000;
options.MaxFuncEvals=1e+5;
options.MaxIters = 1e+5;
options. RelFuncTol=1e-11;
options. StopTol=1e-11;
goodness_X1 = strings(nb_starts,1);
goodness_X = zeros(nb_starts,3-1); %Stores ExitMsg, Fit, F(error for lbfgsb)
for i=1:nb_starts
    if i==1
        
        [Fac_X{i}, ~, out_X{i}] =spca_wopt(X,W,nm_comp,'init','nvecs','alg','ncg','alg_options',options, 'beta', 0.000);
    else
        
        [Fac_X{i}, ~, out_X{i}] =spca_wopt(X,W,nm_comp,'init','random','alg','ncg','alg_options',options,'beta', 0.000);
        
    end
    
    goodness_X1(i) = out_X{i}.ExitFlag;
    goodness_X(i,1) = (1-norm(tensor(times(X,W)-full(Fac_X{i})))^2/norm(tensor(times(X,W)))^2)*100; % fit
    goodness_X(i,2) =norm(tensor(times(X,W)-full(Fac_X{i})))^2; % error
end

%% Get the best Factorization
for i=1:length(goodness_X1)
    good_index(i)=str2num(goodness_X1(i,:));
end
good_flag = find(good_index== 3| good_index==0 );
[er,index]=sort(goodness_X(good_flag,2),'ascend');
Fac = Fac_X{index(1)}; % 
fit=   goodness_X(index(1),1)
%% Get the approximated data
data_T0_appro=Fac.U{1}*diag(Fac.lambda)*(Fac.U{2})'; 
% %% save the approximated data
% save('data_T0_appro.mat','data_T0_appro')

