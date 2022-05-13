
clear all

%%  load dataset
load('Simu_6meta_8time_beta02_set2.mat','X_orig')
% load('Simu_6meta_8time_beta05_set2.mat','X_orig')

%% remove outliers/ subjects with blow-up solution when solving ODE
nr_sub_zeros=find(X_orig.class{1,2}==2); % subjects with blow-up solution
pid_list=str2num(X_orig.label{1});
outlier_index=[nr_sub_zeros]; 
outlier_pid=pid_list(outlier_index);
X_rem=removeoutlier(X_orig,outlier_pid);

%% consider the T0 subtracted data
X=X_rem.data(:,:,2:end)-X_rem.data(:,:,1);
s=size(X);

%% find index for normal and abnormal subjects
sub_normal=find(X_rem.class{1,1}==1);
sub_abnormal=find(X_rem.class{1,1}==2);



%% plot the raw data
labelss = {'Ins_B','GLC_B','Pyr_B','Lac_B','Ala_B','Bhb_B'};
time_value=[0.25 0.5 1 1.5 2 2.5 4];
figure
for i=1:s(2)
    subplot(2,3,i)
    for j=1:length(sub_normal)
        plot(time_value,squeeze(X(sub_normal(j),i,:)),'r')
        hold on
    end
    for j=1:length(sub_abnormal)
        plot(time_value,squeeze(X(sub_abnormal(j),i,:)),'b')
        hold on
    end
    plot(time_value,squeeze(mean(X(sub_normal,i,:),1)),'k','Linewidth',3)
        hold on
    plot(time_value,squeeze(mean(X(sub_abnormal,i,:),1)),'g','Linewidth',3)
    grid on
    xlabel('time(h)')
    ylabel('mmol/L')
    xticks(0:4);
    title(labelss(i))
    set(gca,'FontSize', 15)
end



%% univariate analysis
for k=1:s(3)
for j=1:s(2)
    [h, p_time_meta(k,j), ~, tt]=ttest2(squeeze(X(sub_normal,j,k)), squeeze(X(sub_abnormal,j,k)), 'alpha', 0.05, 'vartype','unequal');
end
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
    X_scal(:,j,:) = temp/rms;    
end
X=tensor(X_scal);


%% plot the preprocessed data
figure
for i=1:s(2)
    subplot(2,3,i)
    for j=1:length(sub_normal)
        plot(time_value,squeeze(X.data(sub_normal(j),i,:)),'r')
        hold on
    end
    for j=1:length(sub_abnormal)
        plot(time_value,squeeze(X.data(sub_abnormal(j),i,:)),'b')
        hold on
    end
    xticks(0:4);
    xlabel('time(h)')
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

%%  Get the best Factorization & fit,TC, CC info 
[er,index]=sort(goodness_X(:,2),'ascend');
Fac = Fac_X{index(1)}; 
out_X_best = out_X{index(1)};
exit_info=goodness_X1(best_F_index(1))
uniqueness=unique_test
fit=out_X_best.Fit
Consistency = corcond(X.data,normalize(Fac,1),[],0)
tc=TC(Fac.U)

%% compute the p-val
for r=1:nm_comp
    r1=Fac.U{1}(sub_normal,r);
    r2=Fac.U{1}(sub_abnormal,r);
    [h, p(r), ~, tt]=ttest2(r1, r2, 'alpha', 0.05, 'vartype','unequal');
end
[pmin,index_pmin]=min(p)


%% plot the factors
Z2={'subjects','metabolites','time'};
xvalue=1:length(Fac.U{1}(:,1));
figure
k=0;
for i=1:3
    
    if i==2
        for j=1:nm_comp
            k=k+1;
            subplot(3,nm_comp,k)
            plot(Fac.U{i}(:,j),'-s','Markersize',6,'LineWidth',2)
            hold on
            ylabel(['comp', num2str(j) ],'FontSize', 10)
            grid on
               ylim([min(min(Fac.U{i})),max(max(Fac.U{i}))])
            set(gca,'xtick',1:length(Fac.U{i}(:,j)),'xticklabel',labelss)
            xtickangle(90)
        end
    elseif i==1
        
        for j=1:nm_comp
            k=k+1;
            subplot(3,nm_comp,k)
            plot(xvalue(sub_normal),Fac.U{i}(sub_normal,j),'rs','MarkerSize',6,...
                'MarkerEdgeColor','r','MarkerFaceColor','r')
            hold on
            plot(xvalue(sub_abnormal),Fac.U{i}(sub_abnormal,j),'bs','MarkerSize',6,...
                'MarkerEdgeColor','b','MarkerFaceColor','b')
            grid on
               ylim([min(min(Fac.U{i})),max(max(Fac.U{i}))])
            ylabel(['comp', num2str(j) ],'FontSize', 10)
            title(Z2{i})
        end
    else
        
        for j=1:nm_comp
            k=k+1;
            subplot(3,nm_comp,k)
            plot(time_value,Fac.U{i}(:,j),'-o','LineWidth',2.4)
            xticks(0:4)
            hold on
            grid on
            ylim([min(min(Fac.U{i})),max(max(Fac.U{i}))])
            ylabel(['comp', num2str(j) ],'FontSize', 10)
           xlabel('time(h)')
        end
    end
    
    
end


%% scatter plot for subjects mode and metabolites mode
% %%
% i=1;
% for j=1:nm_comp
%     for k=j+1:nm_comp
% figure
%       
% 
%             plot(Fac.U{i}(sub_normal,j),Fac.U{i}(sub_normal,k),'rs','MarkerSize',6,...
%                 'MarkerEdgeColor','r','MarkerFaceColor','r')
%             hold on
%             plot(Fac.U{i}(sub_abnormal,j),Fac.U{i}(sub_abnormal,k),'bs','MarkerSize',6,...
%                 'MarkerEdgeColor','b','MarkerFaceColor','b')
%             grid on
%             xlabel(['comp', num2str(j) ],'FontSize', 10)
%             ylabel(['comp', num2str(k) ],'FontSize', 10)
%             set(gca,'Fontsize',13)
%     end
% end
% %%
% i=2;
% for j=1:nm_comp
%     for k=j+1:nm_comp
% figure
% 
% plot(Fac.U{i}(:,j), Fac.U{i}(:,k),'ro','MarkerFaceColor','r');
%         text(Fac.U{i}(:,j), Fac.U{i}(:,k),labelss,'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize', 15,'fontweight', 'bold')
%         xlabel(['PC',num2str(i)])
%         ylabel(['PC',num2str(j)])
%         grid on
%         set(gca,'Fontsize',15)
%               grid on
%             xlabel(['comp', num2str(j) ],'FontSize', 10)
%             ylabel(['comp', num2str(k) ],'FontSize', 10)
%             set(gca,'Fontsize',13)
%     end
% end

