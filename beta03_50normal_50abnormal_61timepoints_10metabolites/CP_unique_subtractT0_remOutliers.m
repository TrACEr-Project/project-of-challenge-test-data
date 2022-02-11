clear all
clc
close all
%% load data

load ('Simu_10meta_61time_beta03.mat','X_orig');
labelss = {'Ins_B','GLC_B','Pyr_B','Lac_B','Ala_B','Glyc_B','Bhb_B','FFA_B',...
           'TG_B','Chol_B'};
%% remove outliers
pid_list=str2num(X_orig.label{1});
outlier_index=[83 85];
outlier_pid=pid_list(outlier_index);
X_rem=removeoutlier(X_orig,outlier_pid);
X=X_rem.data(:,:,2:end)-X_rem.data(:,:,1);
s=size(X);

%% find index for normal and abnormal subjects
 sub_normal=find(X_rem.class{1,1}==1);
sub_abnormal=find(X_rem.class{1,1}==2);

%% plot the raw data
figure
for i=1:s(2)
    subplot(4,3,i)
    for j=sub_normal
            plot(1:s(3),squeeze(X(j,i,:)),'r')
            hold on
    end
    for j=sub_abnormal
            plot(1:s(3),squeeze(X(j,i,:)),'b')
            hold on
    end
    xlim([1 s(3)]);
    title(labelss(i))
    set(gca,'FontSize', 15)    
end

%% preprocessing
%centering across the subjects mode
temp = X(:,:);
temp_centered = temp - repmat(mean(temp),size(temp,1),1);
X_centered = reshape(temp_centered, size(X));
%scaling within the second(metabolites) mode - using root mean square
for j=1:size(X_centered,2)
    temp = squeeze(X_centered(:,j,:));
    rms = sqrt(mean((temp(:).^2)));
    X_scal(:,j,:) = temp/rms;           %scaling by rms
end
X=tensor(X_scal);


%% plot the preprocessed data
figure
for i=1:s(2)
    subplot(4,3,i)
    for j=sub_normal
            plot(1:s(3),squeeze(X.data(j,i,:)),'r')
            hold on
    end
    for j=sub_abnormal
            plot(1:s(3),squeeze(X.data(j,i,:)),'b')
            hold on
    end
    xlim([1 s(3)]);
    title(labelss(i))
    set(gca,'FontSize', 15)    
end



%% CP MODEL
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
goodness_X = zeros(nb_starts,2); %Stores ExitMsg, Fit, F(error for lbfgsb)
% Fac_X = cell(nb_starts,1);
% out_X = cell(nb_starts,1);
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
%% uniqueness test
% 0 -> NOT unique
% 1 -> Unique
% 2 -> Inconclusive, need more random starts
good_flag = find(goodness_X1(:) == 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL.' | goodness_X1(:) == 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH.');
if length(good_flag)>=1
    F_round = round(goodness_X(good_flag,2),7);
    best_F_index = good_flag(F_round == min(F_round));
    if length(best_F_index) < 2 
        F_round = round(goodness_X(good_flag,2),5); 
        best_F_index = good_flag(F_round == min(F_round));
    end
else
    F_round = round(goodness_X(:,2),10);
    best_F_index = find(F_round == min(F_round));
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

%%  Get the best Factorization
[er,index]=sort(goodness_X(:,2),'ascend');
Fac = Fac_X{index(1)}; % best Fac
out_X_best = out_X{index(1)};
exit_info=goodness_X1(best_F_index(1))
uniqueness=unique_test
fit=out_X_best.Fit
Consistency = corcond(X.data,normalize(Fac,1),[],0)
tc=TC(Fac.U)
%% compute the similarity score in each mode
c=zeros(nm_comp,nm_comp); % subjects mode
S=zeros(nm_comp,nm_comp); % metabolites mode
t=zeros(nm_comp,nm_comp); % time mode
for i=1:nm_comp
    for j=1:nm_comp
        a1=Fac.U{1}(:,i);a2=Fac.U{1}(:,j);
        c(i,j)=a1'*a2/norm(a1)/norm(a2);
        b1=Fac.U{2}(:,i);b2=Fac.U{2}(:,j);
        S(i,j)=b1'*b2/norm(b1)/norm(b2);
        c1=Fac.U{3}(:,i);c2=Fac.U{3}(:,j);
        t(i,j)=c1'*c2/norm(c1)/norm(c2);
    end
end

%% compute the p-val

for r=1:nm_comp
    r1=Fac.U{1}(sub_normal,r);
    r2=Fac.U{1}(sub_abnormal,r);
        [h, p(r), ~, tt]=ttest2(r1, r2, 'alpha', 0.05, 'vartype','unequal');
end
[pmin,index_pmin]=min(p)


%%  plot
s=size(X);
Z2={'subjects','metabolites','time'};
Leglab = {'Comp1', 'Comp2','Comp3', 'Comp4','Comp5', 'Comp6','Comp7','Comp8', 'Comp9','Comp10'};

%%% plot the factors
figure 
k=0;
for i=1:3
    
    if i==2
        for j=1:nm_comp
            k=k+1;
            subplot(3,nm_comp,k)
            plot(Fac.U{i}(:,j),'-o','LineWidth',2.4)
            hold on
            ylabel(['comp', num2str(j) ],'FontSize', 10)
            grid on
            set(gca,'xtick',1:length(Fac.U{i}(:,j)),'xticklabel',labelss)
            xtickangle(90)   
        end
    else
        for j=1:nm_comp
            k=k+1;
            subplot(3,nm_comp,k)
            plot(Fac.U{i}(:,j),'-o','LineWidth',2.4)
            hold on
            grid on
            ylabel(['comp', num2str(j) ],'FontSize', 10)
            title(Z2{i})
        end
    end

    
end


% %% check outlier
% check_outliers_CP(X, Fac)