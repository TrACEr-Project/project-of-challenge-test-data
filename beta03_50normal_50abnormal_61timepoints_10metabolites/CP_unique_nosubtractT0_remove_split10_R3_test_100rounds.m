
clear all
clc
% addpath('/Users/ll/Documents/MATLAB/toolbox/PLS_Toolbox_891')
%% load data

load ('Simu_10meta_61time_beta03.mat','X_orig');
labelss = {'Ins_B','GLC_B','Pyr_B','Lac_B','Ala_B','Glyc_B','Bhb_B','FFA_B',...
    'TG_B','Chol_B'};
%% remove outliers
pid_list=str2num(X_orig.label{1});
outlier_index=[83 85];
outlier_pid=pid_list(outlier_index);
Y_rem=removeoutlier(X_orig,outlier_pid);

%% find index for normal and abnormal subjects
sub_normal=find(Y_rem.class{1,1}==1);
sub_abnormal=find(Y_rem.class{1,1}==2);
index_perm=[sub_normal,sub_abnormal]; % make sure that the normal (appear first) are together and abnormal are together
%% kk (10) round of split10
for kk=1:10
    S_perm=[randperm(length(sub_normal)),randperm(length(sub_abnormal))+length(sub_normal)];
    Y_split_left=cell(1,10);
    Results=cell(1,10);
    for i=1:10
        S_perm_rem=S_perm(i:10:end);
        pid_list=str2num(Y_rem.label{1});
        remove_pid=pid_list(S_perm_rem);
        Y_split_left{i}=removeoutlier(Y_rem,remove_pid')
        Results{i}.Sub_rem=S_perm_rem;  
    end
    for ii=1:length(Y_split_left)
        X=tensor(Y_split_left{ii}.data);
        % preprocessing
        %centering across the condition mode
        XX=X.data;
        temp = XX(:,:);
        temp_centered = temp - repmat(mean(temp),size(temp,1),1);
        XXX = reshape(temp_centered, size(XX));
        %scaling in the second mode - using std or root mean square
        X=tensor(XXX);
        for j=1:size(X,2)
            temp = squeeze(X.data(:,j,:));
            rms = sqrt(mean((temp(:).^2)));
            XX(:,j,:) = temp/rms;           %scaling by rms
        end
        Xpre=tensor(XX);
        Results{ii}.Xpre=Xpre;
        % CP model
        X=tensor(Xpre);
        nb_starts =20;
        nm_comp=3;
        optionsCP.factr=1e6;
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
                
                [Fac_X{i}, ~, out_X{i}] =cp_opt(X,nm_comp,'init','nvecs','lower',Low,'opt_option',optionsCP);
            else
                
                [Fac_X{i}, ~, out_X{i}] =cp_opt(X,nm_comp,'init','randn','lower',Low,'opt_option',optionsCP);
                
            end
            
            goodness_X1(i) = out_X{i}.ExitMsg;
            goodness_X(i,1) = out_X{i}.Fit;
            goodness_X(i,2) = out_X{i}.OptOut.err(end,1);
        end
        Results{ii}.Fac_X=Fac_X;
        Results{ii}.cov_info=goodness_X1;
        Results{ii}.er=goodness_X(:,2);
        
        %%uniqueness test
        % 0 -> NOT unique
        % 1 -> Unique
        % 2 -> Inconclusive, need more random starts
        good_flag = find(goodness_X1(:) == 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL.' | goodness_X1(:) == 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH.');
        % good_flag =1:nb_starts;
        if length(good_flag)>=1
            F_round = round(goodness_X(good_flag,2),6);
            best_F_index = good_flag(F_round == min(F_round));
            if length(best_F_index) < 2
                F_round = round(goodness_X(good_flag,2),2);
                best_F_index = good_flag(F_round == min(F_round));%Finds best F values
            end
            %     flag = 0;
        else
            F_round = round(goodness_X(:,2),10);
            best_F_index = find(F_round == min(F_round));
            %     flag = 1;
        end
        
        eps = .05;
        if length(best_F_index)==1
            unique_test = 2;
            disp('Need more random starts to determine uniqueness')
            worst_check = 0;
        elseif length(best_F_index) > 1
            check_matrix = zeros(length(best_F_index));
            for i = 1:length(best_F_index)
                for j = 1:length(best_F_index)
                    check_matrix(i,j) = score(Fac_X{best_F_index(j)},Fac_X{best_F_index(i)},'lambda_penalty',false);
                end
            end
            worst_check = min(min(check_matrix));
            if worst_check < (1-eps) %Checks to see if factors are the same if F is
                unique_test = 0;
            else
                unique_test = 1;
            end
        end
        unique_test
        nr_best_F_index=length(best_F_index);
        Results{ii}.uniqueness=unique_test;
        Results{ii}.nr_best_F_index=nr_best_F_index;
        
        %%Get the best Factorization
        [er,index]=sort(goodness_X(:,2),'ascend');
        Fac = Fac_X{index(1)}; % best Fac
        out_X_best = out_X{index(1)};
        exit_info=goodness_X1(best_F_index(1))
        uniqueness=unique_test
        fit=out_X_best.Fit
        Consistency = corcond(X.data,normalize(Fac,1),[],0)
        tc=TC(Fac.U)
        %%record these results
        Results{ii}.Fac=Fac;
        Results{ii}.fit=fit;
        Results{ii}.tc=tc;
        Results{ii}.CC=Consistency;
        %%compute the p_val
        for r=1:nm_comp
            index1=find(Y_split_left{ii}.class{1,1}==1);
            index2=find(Y_split_left{ii}.class{1,1}==2);
            r1=Fac.U{1}(index2,r);
            r2=Fac.U{1}(index1,r);
            [h, pval(r), ~, tt]=ttest2(r1, r2, 'alpha', 0.05, 'vartype','unequal');
        end
        pmin=min(pval)
        Results{ii}.pval=pmin;
    end
    
    % pick out the metabolites mode and time mode 
    for ii=1:10
        Fac_rem_submode{ii}=ktensor(Results{ii}.Fac.lambda,Results{ii}.Fac.U{2},Results{ii}.Fac.U{3});
    end
    
    for ii=1:10
        [scores_M_T(ii), Permuted_Fac_M_T{ii}] = score(Fac_rem_submode{ii},Fac_rem_submode{1},'lambda_penalty',false);
    end
    R3_Results_CP_all_centered_test_100round.scoresMT{kk}=scores_M_T;
    R3_Results_CP_all_centered_test_100round.worstscoresMT(kk)=min(scores_M_T);
    R3_Results_CP_all_centered_test_100round.results{kk}=Results;
end
%%
save ('R3_Results_CP_all_centered_test_100round.mat','R3_Results_CP_all_centered_test_100round')


%% plot factors
load('R3_Results_CP_all_centered_test_100round.mat','R3_Results_CP_all_centered_test_100round')

for kk=1:10
    Results=R3_Results_CP_all_centered_test_100round.results{kk};
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
title('3-comp CP model')
set(gca,'Fontsize',18)

%  Permuted_Fac_M_T{1}.U{1}(:,3)=- Permuted_Fac_M_T{1}.U{1}(:,3);
nm_comp=3;
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
labelss = {'Ins_B','GLC_B','Pyr_B','Lac_B','Ala_B','Glyc_B','Bhb_B','FFA_B',...
    'TG_B','Chol_B'};

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
            ylim([-0.5 1])
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





