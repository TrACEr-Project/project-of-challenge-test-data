% This script demonstrates the replicability check of the CP model to the real data set with only male subjects and 162 measurements 

% We use Tensor Toolbox as well as the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C
% In addition, parts of the scripts may require the dataset object (https://eigenvector.com/software/dataset-object/), publically available.
% For core consistency computation, we also use the corcord function from the Nway toolbox (http://www.models.life.ku.dk/nwaytoolbox).


clear all
clc
close all

%% add auxilary functions to path
addpath(genpath('.\functions'))
%% add other apckages to your path!
addpath(genpath('...\tensor_toolbox-v3.1')) %Tensor toolbox is needed;  MATLAB Tensor Toolbox. Copyright 2017, Sandia Corporation, http://www.tensortoolbox.org/
addpath(genpath('...\L-BFGS-B-C-master')) % LBFGS-B implementation is needed; download here: https://github.com/stephenbeckr/L-BFGS-B-C
addpath(genpath('...\nway331')) % Nway toolbox is needed for computing core consistency; download here: http://www.models.life.ku.dk/nwaytoolbox
addpath(genpath('...\dataset')) % dataset object is needed; download here: https://eigenvector.com/software/dataset-object/

%%
load('NMR+Gene+SNP+Metav3_July4_2022.mat') %% load the real data

%% remove outliers
pid_list = str2num(NMR.label{1});
outlier_pid =[142,79,342,250,90,335,312];
NMR_remove=removesubject(NMR,outlier_pid);

%% remove the metabolites with lots of missing values
NMR_remove=removeisnan(NMR_remove);

%% select the large set of metabolites (162 measurements)
selected_metabolite_index = find(NMR_remove.class{3,1}==1 | NMR_remove.class{3,1}==2);
NMR_remove= NMR_remove(:,:,selected_metabolite_index);

%% choose male(2)/female(1)
selected_gender_index =find(NMR_remove.class{1,1}==2);
NMR_remove= NMR_remove(selected_gender_index,:,:);


%% BMI index
index_Under=find(NMR_remove.class{1,2}==4);
index_normal=find(NMR_remove.class{1,2}==1);
index_obesity=find(NMR_remove.class{1,2}==2);
index_over=find(NMR_remove.class{1,2}==3);
sub_normal=[index_Under,index_normal]; % Low BMI group
sub_abnormal=[index_over,index_obesity]; % High BMI group

% rearrange the subjects so that the Low BMI subjects appears first followed by High BMI subjects
index_perm=[sub_normal, sub_abnormal]; 
% data considered, before split
Y_rem=NMR_remove(index_perm,:,:);

% kk (=10) round of split10
for kk=1:10
    S_perm=[randperm(length(sub_normal)),randperm(length(sub_abnormal))+length(sub_normal)];
    Y_split_left=cell(1,10);
    Results=cell(1,10);
    for i=1:10
        S_perm_rem=S_perm(i:10:end);
        pid_list=str2num(Y_rem.label{1});
        remove_pid=pid_list(S_perm_rem);
        Y_split_left{i}=removesubject(Y_rem,remove_pid');
        Results{i}.Sub_rem=S_perm_rem;
    end
    for ii=1:length(Y_split_left)
        X=tensor(Y_split_left{ii}.data(:,2:end,:)-Y_split_left{ii}.data(:,1,:));
        % preprocess the data
        s=size(X);
        X_center=X.data-repmat(nanmean(X.data,1),s(1),1);
        for iii=1:s(3)
            temp=X_center(:,:,iii);
            X_pre(:,:,iii)=temp/sqrt(nanmean(temp.^2,'all'));
        end
        % 
        Xpre=permute(X_pre, [1 3 2]);
        Results{ii}.Xpre=Xpre;
        X=tensor(Xpre);
        % CP model (CP_wopt) 
        W=ones(size(X));
        W(find(isnan(X.data)))=0;
        W=tensor(W);
        X(find(isnan(X.data)))=0;
        nb_starts =20;
        nm_comp=2;
        optionsCP.factr=1e3;
        optionsCP.maxIts = 10000;
        optionsCP.maxTotalITs=50000;
        optionsCP.printEvery  = 10000;
        Low{1}=-Inf*ones(size(X,1),nm_comp);
        Low{2}=-Inf*ones(size(X,2),nm_comp);
        Low{3}=-Inf*ones(size(X,3),nm_comp);
        goodness_X1 = strings(nb_starts,1);
        goodness_X = zeros(nb_starts,3-1); %Stores ExitMsg, Fit, F(error for lbfgsb)
        %Call cp_opt to fit the CP model
        for i=1:nb_starts
            if i==1
                
                [Fac_X{i}, ~, out_X{i}] =cp_wopt(X,W,nm_comp,'init','nvecs','lower',Low,'opt_option',optionsCP);
            else
                
                [Fac_X{i}, ~, out_X{i}] =cp_wopt(X,W,nm_comp,'init','randn','lower',Low,'opt_option',optionsCP);
                
            end
            
            goodness_X1(i) = out_X{i}.ExitMsg;
            goodness_X(i,1) = (1-norm(tensor(times(X,W)-full(Fac_X{i})))^2/norm(tensor(times(X,W)))^2)*100; % fit
            goodness_X(i,2) =norm(tensor(times(X,W)-full(Fac_X{i})))^2; % error
        end
        Results{ii}.Fac_X=Fac_X;
        Results{ii}.cov_info=goodness_X1;
        Results{ii}.er=goodness_X(:,2);
        
        %Get the best Factorization
        [er,index]=sort(goodness_X(:,2),'ascend');
        Fac = Fac_X{index(1)}; % best Fac
        exit_info=goodness_X1(index(1));
        fit=goodness_X(index(1),1);
        Consistency = corcond(X.data,normalize(Fac,1),[],0);
        tc=TC(Fac.U);
        %record these results
        Results{ii}.Fac=Fac;
        Results{ii}.fit=fit;
        Results{ii}.tc=tc;
        Results{ii}.CC=Consistency;
        
    end
    Results_all{kk}.Results=Results;

end



%% 
%   compute the factor match score (FMS) by the two sets of factors (only use the factors on the metabolites and time modes) 
%   extracted from any two datasets from all splits 

for kk=1:10
    for ii=1:10
        Fac_rem_submode{(kk-1)*10+ii}=ktensor(Results_all{kk}.Results{ii}.Fac.lambda,...
                                              Results_all{kk}.Results{ii}.Fac.U{2},...
                                              Results_all{kk}.Results{ii}.Fac.U{3});
    end
end
m=0;
for k=1:length(Fac_rem_submode)
    for l=k+1:length(Fac_rem_submode)
        m=m+1;
        [scores_M_T(m), Permuted_Fac_M_T{m}] = score(Fac_rem_submode{k},Fac_rem_submode{l},'lambda_penalty',false);
    end
end

f=figure
h=histogram(scores_M_T,[0:0.1:1],'Normalization','probability');
xlim([0 1])
xticks([0 0.5 0.9 1.0])
xlabel('FMS')
ylabel('Frequency')
set(gca,'Fontsize',18)
ylim([0 1.0])
yticks([0:0.1:1.0])
set(gca,'Fontsize',18)



