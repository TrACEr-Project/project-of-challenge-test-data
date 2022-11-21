% This script shows an example to check the replicability of the CP model to the simulated datasets

% We use Tensor Toolbox as well as the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C
% In addition, parts of the scripts may require the dataset object (https://eigenvector.com/software/dataset-object/), publically available.
% For core consistency computation, we also use the corcord function from the Nway toolbox (http://www.models.life.ku.dk/nwaytoolbox).

clear all

%%  load dataset
load('Simu_6meta_8time_alpha02_IRM_balance.mat','X_orig')
% load('Simu_6meta_8time_alpha04_IRM_balance.mat','X_orig')
% load('Simu_6meta_8time_alpha02_IRM_unbalance.mat','X_orig')
% load('Simu_6meta_8time_alpha04_IRM_unbalance.mat','X_orig')
% load('Simu_6meta_8time_alpha02_betacell_balance.mat','X_orig')
% load('Simu_6meta_8time_alpha04_betacell_balance.mat','X_orig')
% load('Simu_6meta_8time_alpha02_betacell_unbalance.mat','X_orig')
% load('Simu_6meta_8time_alpha04_betacell_unbalance.mat','X_orig')

%% remove subjects with blow-up solution when solving ODE / remove outliers
nr_sub_zeros=find(X_orig.class{1,2}==2); % subjects with blow-up solution
pid_list=str2num(X_orig.label{1});
outlier_index=[nr_sub_zeros]; 
outlier_pid=pid_list(outlier_index);
Y_rem=removesubject(X_orig,outlier_pid);

%% rearrange the data so that all the normal subjects appear first 
sub_normal=find(Y_rem.class{1,1}==1);
sub_abnormal=find(Y_rem.class{1,1}==2);
index_perm=[sub_normal,sub_abnormal]; 

%% data before split
Y_rem=Y_rem(index_perm,:,:);

%% in total, we do split10 10 rounds; kk-- the kk_th round;
for kk=1:10
    
    S_perm=[randperm(length(sub_normal)),randperm(length(sub_abnormal))+length(sub_normal)];
    Y_split_left=cell(1,10); % to store the 10 splits in each split10
    Results=cell(1,10);% to store the results in each round of split10
    
    
    % randomly remove 1/10 of subjects
    for i=1:10
        S_perm_rem=S_perm(i:10:end);
        pid_list=str2num(Y_rem.label{1});
        remove_pid=pid_list(S_perm_rem);
        Y_split_left{i}=removesubject(Y_rem,remove_pid') % the data left after removing 1/10 part
        Results{i}.Sub_rem=S_perm_rem;  % record the position of the removed data
    end
    
    
    % Preprocess the data and run CP model for each split in split10
    for ii=1:length(Y_split_left)
        X=tensor(Y_split_left{ii}.data);
        % preprocessing
        %centering across the condition mode
        XX=X.data;
        temp = XX(:,:);
        temp_centered = temp - repmat(mean(temp),size(temp,1),1);
        XXX = reshape(temp_centered, size(XX));
        %scaling in the second mode - using root mean square
        X=tensor(XXX);
        for j=1:size(X,2)
            temp = squeeze(X.data(:,j,:));
            rms = sqrt(mean((temp(:).^2)));
            XX(:,j,:) = temp/rms; %scaling by rms
        end
        Xpre=tensor(XX);
        Results{ii}.Xpre=Xpre; % record the preprocessed data in each split
        %run CP model
        X=tensor(Xpre);
        nb_starts =20;
        nm_comp=4;
        optionsCP.factr=1e6;
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
                
                [Fac_X{i}, ~, out_X{i}] =cp_opt(X,nm_comp,'init','nvecs','lower',Low,'opt_option',optionsCP);
            else
                
                [Fac_X{i}, ~, out_X{i}] =cp_opt(X,nm_comp,'init','randn','lower',Low,'opt_option',optionsCP);
                
            end
            
            goodness_X1(i) = out_X{i}.ExitMsg;
            goodness_X(i,1) = out_X{i}.Fit;
            goodness_X(i,2) = out_X{i}.OptOut.err(end,1);
        end
        % record the results from all initializations for future check
        Results{ii}.Fac_X=Fac_X;
        Results{ii}.cov_info=goodness_X1;
        Results{ii}.er=goodness_X(:,2);
        unique_test = unique_test_CP(Fac_X, goodness_X(:,2), goodness_X1)
        Results{ii}.uniqueness=unique_test;
        
        %Get info for the best Factorization
        [er,index]=sort(goodness_X(:,2),'ascend');
        Fac = Fac_X{index(1)}; % best Fac
        out_X_best = out_X{index(1)};
        exit_info=goodness_X1(index(1))
        uniqueness=unique_test
        fit=out_X_best.Fit
        Consistency = corcond(X.data,normalize(Fac,1),[],0)
        tc=TC(Fac.U)
        for r=1:nm_comp %compute the p_val for each comp in the best Fac
            index1=find(Y_split_left{ii}.class{1,1}==1);
            index2=find(Y_split_left{ii}.class{1,1}==2);
            r1=Fac.U{1}(index2,r);
            r2=Fac.U{1}(index1,r);
            [h, pval(r), ~, tt]=ttest2(r1, r2, 'alpha', 0.05, 'vartype','unequal');
        end
        pmin=min(pval)
        
        
        %%record these results
        Results{ii}.Fac=Fac;
        Results{ii}.fit=fit;
        Results{ii}.tc=tc;
        Results{ii}.CC=Consistency;
        Results{ii}.pval=pmin;
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

f=figure;
h=histogram(scores_M_T,[0:0.1:1],'Normalization','probability');
xlim([0 1])
xticks([0 0.5 0.9 1.0])
xlabel('FMS')
ylabel('Frequency')
set(gca,'Fontsize',18)
ylim([0 1.0])
yticks([0:0.1:1.0])
set(gca,'Fontsize',18)






