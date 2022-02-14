clear all
clc
close all
% need to add tensor toolbox
% addpath('/Users/ll/Documents/MATLAB/toolbox/tensor_toolbox_v3.1')
% addpath('/Users/ll/Documents/MATLAB/toolbox/PLS_Toolbox_891') % for datset

%% load the real data
load('NMR+Gene+SNP+Metav3_Feb11.mat')

%% pick out the common metabolites with the simulated data
% The intersected metabolites are {'Ins','GLC','Pyr','Lac','Ala','Bhb'};
NMR=select_metabolites_simul(NMR); % If you use other metabolites,just remove this line

%% remove subjects and metabolites that have too many nans
NMR=removeisnan(NMR);

%% reorder the metabolites as {'Ins','GLC','Pyr','Lac','Ala','Bhb'};
index=[6 2 4 3 1 5 ];
NMR=NMR(:,:,index);
labelss = {'Ins','GLC','Pyr','Lac','Ala','Bhb'};

%% parameters
intersec_gene=1;%% 0 means not intersected data;1 means intersected data
gender='male'; %% options are 'female','male','full'
meta_type=0; %% options are 1,2 and others;1 means Age_selected(small) sets;2 means Gozde select(large) sets.
data_type='full dynamic'; %% options are: 'T0', 'T0 corrected', 'full dynamic'

%% select data (gene-intersection,male/female/full) from the real data
NMR=select_data(NMR,Gene,intersec_gene,gender,meta_type);

%% remove outliers
pid_list = str2num(NMR.label{1});
index = [];
index=[ 36 67 70 110 122];
outlier_pid =  [pid_list(index)];
NMR=removeoutlier(NMR,outlier_pid);

%% find index for normal and abnormal subjects
sub_normal=find(NMR.class{1,2}==4|NMR.class{1,2}==1);
sub_abnormal=find(NMR.class{1,2}==2|NMR.class{1,2}==3);
index_perm=[sub_normal,sub_abnormal]; % make sure that the normal (appear first) are together and abnormal are together
%% data considered, before split
Y_rem=NMR(index_perm,:,:);


%% kk (10) round of split10
for kk=1:10
    S_perm=[randperm(length(sub_normal)),randperm(length(sub_abnormal))+length(sub_normal)];
    Y_split_left=cell(1,10);
    Results=cell(1,10);
    for i=1:10
        S_perm_rem=S_perm(i:10:end);
        pid_list=str2num(Y_rem.label{1});
        remove_pid=pid_list(S_perm_rem);
        Y_split_left{i}=removeoutlier(Y_rem,remove_pid)
        Results{i}.Sub_rem=S_perm_rem;
    end
    for ii=1:length(Y_split_left)
        X=Y_split_left{ii}.data;
        % preprocessing
        X=tensor(preprocess(X));
        Results{ii}.Xpre=X;
        %%cp_wopt model
        % create index matrix
        W=ones(size(X));
        W(find(isnan(X.data)))=0;
        W=tensor(W);
        X(find(isnan(X.data)))=0;
        % cp_wopt
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
        Fac = Fac_X{index(1)};
        uniqueness=unique_test
        fit=   goodness_X(index(1),1)
        exit_info=goodness_X1(best_F_index(1))
        Consistency = corcond(X.data,normalize(Fac,1),[],0)
        tc=TC(Fac.U)
        
        %%record these results
        Results{ii}.Fac=Fac;
        Results{ii}.fit=fit;
        Results{ii}.tc=tc;
        Results{ii}.CC=Consistency;
        %%compute the p_val
        
        for r=1:nm_comp
            index1=find(Y_split_left{ii}.class{1,2}==4|Y_split_left{ii}.class{1,2}==1);
            index2=find(Y_split_left{ii}.class{1,2}==2|Y_split_left{ii}.class{1,2}==3);
            r1=Fac.U{1}(index2,r);
            r2=Fac.U{1}(index1,r);
            [h, pval(r), ~, tt]=ttest2(r1, r2, 'alpha', 0.05, 'vartype','unequal');
        end
        Results{ii}.pval=min(pval);
        Results{ii}.nrabnormal=length(index2);
        Results{ii}.nrnormal=length(index1);
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

%%
for i=1:10
    nr_normal(i)=Results{i}.nrnormal;
    nr_abnormal(i)=Results{i}.nrabnormal;
end


