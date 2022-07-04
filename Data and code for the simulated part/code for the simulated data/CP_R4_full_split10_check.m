

clear all

%%  load dataset
load('Simu_6meta_8time_beta02_IRM_balance.mat','X_orig')
% load('Simu_6meta_8time_beta04_IRM_balance.mat','X_orig')
% load('Simu_6meta_8time_beta02_IRM_unbalance.mat','X_orig')
% load('Simu_6meta_8time_beta04_IRM_unbalance.mat','X_orig')
% load('Simu_6meta_8time_beta02_betacell_balance.mat','X_orig')
% load('Simu_6meta_8time_beta04_betacell_balance.mat','X_orig')
% load('Simu_6meta_8time_beta02_betacell_unbalance.mat','X_orig')
% load('Simu_6meta_8time_beta04_betacell_unbalance.mat','X_orig')

%% remove subjects with blow-up solution when solving ODE / remove outliers
nr_sub_zeros=find(X_orig.class{1,2}==2); % subjects with blow-up solution
pid_list=str2num(X_orig.label{1});
outlier_index=[nr_sub_zeros]; 
outlier_pid=pid_list(outlier_index);
Y_rem=removeoutlier(X_orig,outlier_pid);

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
    Results=cell(1,10);% to store the results in each split10
    
    
    % randomly remove 1/10 of subjects
    for i=1:10
        S_perm_rem=S_perm(i:10:end);
        pid_list=str2num(Y_rem.label{1});
        remove_pid=pid_list(S_perm_rem);
        Y_split_left{i}=removeoutlier(Y_rem,remove_pid') % the data left after removing 1/10 part
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
        % record the results from all initialization for future check
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
    
    R4_Results_CP_all_10split10.results{kk}=Results; % to stop the results in each round of split10
end


% %% 
% save ('R4_Results_CP_all_10split10.mat','R4_Results_CP_all_10split10')
% 
% 
% %% plot factors
% load('R4_Results_CP_all_10split10.mat','R4_Results_CP_all_10split10')

for kk=1:length(R4_Results_CP_all_10split10)  % loop through each round of split10
    Results=R4_Results_CP_all_10split10.results{kk};
    % pick out the metabolites mode and time mode
    for ii=1:10
        Fac_rem_submode{(kk-1)*10+ii}=ktensor(Results{ii}.Fac.lambda,Results{ii}.Fac.U{2},Results{ii}.Fac.U{3});
    end
end
for k=1:length(Fac_rem_submode)
    
    [scores_M_T(k), Permuted_Fac_M_T{k}] = score(Fac_rem_submode{k},Fac_rem_submode{1},'lambda_penalty',false);
    
end

figure
h=histogram(scores_M_T,'Normalization','probability');
xlim([0 1])
xlabel('Factor match score')
ylabel('Frequency')
title('4-comp CP model')
set(gca,'Fontsize',18)

% flip the sign to make the factors comparable
nm_comp=4;
for k=1:length(Fac_rem_submode)
    for j=1:nm_comp
        if Permuted_Fac_M_T{k}.U{2}(:,j)'*Permuted_Fac_M_T{1}.U{2}(:,j)<0
            Permuted_Fac_M_T{k}.U{2}(:,j)=-Permuted_Fac_M_T{k}.U{2}(:,j);
        end
        if Permuted_Fac_M_T{k}.U{1}(:,j)'*Permuted_Fac_M_T{1}.U{1}(:,j)<0
            Permuted_Fac_M_T{k}.U{1}(:,j)=-Permuted_Fac_M_T{k}.U{1}(:,j);
        end
    end
    
    
end


figure
labelss = {'Ins_B','GLC_B','Pyr_B','Lac_B','Ala_B','Bhb_B'};
kk=0;
for i=1:2
    for j=1:nm_comp
        kk=kk+1;
        if i==1
            subplot(2,nm_comp,kk)
            for k=1:length(Permuted_Fac_M_T)
                plot(Permuted_Fac_M_T{k}.U{i}(:,j))
                hold on
            end
            ylim([-1 1])
            set(gca,'xtick',1:length(Permuted_Fac_M_T{k}.U{i}(:,j)),'xticklabel',labelss)
            xtickangle(90)
            ylabel(['comp', num2str(j)])
            grid on
        else
            subplot(2,nm_comp,kk)
            for k=1:length(Permuted_Fac_M_T)
                plot(Permuted_Fac_M_T{k}.U{i}(:,j))
                hold on
            end
            ylabel(['comp', num2str(j)])
            grid on
            
        end
    end
end





