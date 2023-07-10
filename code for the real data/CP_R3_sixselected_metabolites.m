% This script shows how to fit a CP model to the real data set with only male subjects and the selected 162 measurements 

% We use Tensor Toolbox as well as the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C
% In addition, parts of the scripts may require the dataset object (https://eigenvector.com/software/dataset-object/), publically available.
% For core consistency computation, we also use the corcord function from the Nway toolbox (http://www.models.life.ku.dk/nwaytoolbox).


clear all
clc
close all
%% add auxilary functions to path
addpath(genpath('./functions'))
%% add other apckages to path
addpath(genpath('.../tensor_toolbox-v3.1')) %Tensor toolbox is needed;  MATLAB Tensor Toolbox. Copyright 2017, Sandia Corporation, http://www.tensortoolbox.org/
addpath(genpath('.../L-BFGS-B-C-master')) % LBFGS-B implementation is needed; download here: https://github.com/stephenbeckr/L-BFGS-B-C
addpath(genpath('.../nway331')) % Nway toolbox is needed for computing core consistency; download here: http://www.models.life.ku.dk/nwaytoolbox
addpath(genpath('.../dataset')) % dataset object is needed; download here: https://eigenvector.com/software/dataset-object/

%%
load('data.mat') %% load the real NMR data (dataset format)

%% remove outliers
outlier_pid =[]; % remove outliers
NMR_remove=removesubject(NMR,outlier_pid);

%% remove the metabolites with lots of missing values
NMR_remove=removeisnan(NMR_remove);

%% select the six metabolites showing up in the metabolic model
selected_metabolite_index = [249 70 72 71 61 74];
NMR_remove= NMR_remove(:,:,selected_metabolite_index);


%% Get the raw data
X_slect = NMR_remove.data;

%% preprocessing (centering accross subject mode and scaling within metabolite mode)
s=size(X_slect);
X_center=X_slect-repmat(nanmean(X_slect,1),s(1),1);
for i=1:s(3)
    temp=X_center(:,:,i);
    X_pre(:,:,i)=temp/sqrt(nanmean(temp.^2,'all'));
end
X=X_pre;
X=permute(X, [1 3 2]);
X=tensor(X);

%%  cp_wopt model
% create missing-value index matrix
W=ones(size(X));
W(find(isnan(X.data)))=0;
W=tensor(W);
X(find(isnan(X.data)))=0;
% set parameters for cp_wopt
nb_starts =32; % number of random initializations
nm_comp=3;
optionsCP.factr=1e3;
optionsCP.maxIts = 10000;
optionsCP.maxTotalITs=50000;
optionsCP.printEvery  = 10000;
Low{1}=-Inf*ones(size(X,1),nm_comp);
Low{2}=-Inf*ones(size(X,2),nm_comp);
Low{3}=-Inf*ones(size(X,3),nm_comp);
goodness_X1 = strings(nb_starts,1);
goodness_X = zeros(nb_starts,3-1); %Stores ExitMsg, Fit, F(error for lbfgsb)
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


%% numerically check the uniqueness of the CP factorization
unique_test=unique_test_CP(Fac_X,goodness_X(:,2),goodness_X1);
good_flag = find(goodness_X1(:) == 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL.' | goodness_X1(:) == 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH.');

%% Get the best CP factorization
[er,index]=sort(goodness_X(good_flag,2),'ascend');
Fac = Fac_X{index(1)};

%% Check fit, exit_info, core consistency, tucker congruency for the CP decomposition
fit= goodness_X(index(1),1);
exit_info=goodness_X1(index(1));
Consistency = corcond(X.data,normalize(Fac,1),[],0);
tc=TC(Fac.U);



%% plot of the CP factors
labelss = {'Ins','Glc','Pyr','Lac','Ala','Bhb'};
f=figure;
k=0;
for ii=1:3
    for jj = 1:nm_comp
        k=k+1;
        subplot(3, nm_comp, k)
        x_value = 1:length(Fac.U{ii}(:,1));
        y_value =Fac.U{ii}(:,jj);
        if ii==2
            plot(x_value,y_value,'-s','MarkerSize',6);
            ylabel(['b_',num2str(jj)],'fontweight','bold')
            ylim([min(min(Fac.U{ii})),max(max(Fac.U{ii}))])
            xlim([min(x_value) max(x_value)])
            set(gca,'xtick',x_value,'xticklabel',labelss)
            set(gca,'Fontsize',18)
        elseif ii==1
            
            plot(x_value,y_value,'s')
            ylabel(['a_',num2str(jj)],'fontweight','bold')
            xlabel(['subjects'])
            ylim([min(min(Fac.U{ii})),max(max(Fac.U{ii}))])
            xlim([min(x_value) max(x_value)])
            set(gca,'Fontsize',18)
        else
            time_values=[0 0.25 0.5 1 1.5 2 2.5 4];
            plot(time_values,y_value,'-','Linewidth',3)
            xticks(0:4)
            ylabel(['c_',num2str(jj)],'fontweight','bold')
            xlabel(['time(h)'])
            ylim([min(min(Fac.U{ii})),max(max(Fac.U{ii}))])
            set(gca,'Fontsize',18)
        end
        
    end
end




