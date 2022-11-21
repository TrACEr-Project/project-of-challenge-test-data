% This script shows how to fit the PCA model using the svd algorithm to the simulated data set 

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
X_rem=removesubject(X_orig,outlier_pid);


%% consider the T0 data
T0=squeeze(X_rem.data(:,:,1));

%% find index for normal and abnormal subjects
sub_normal=find(X_rem.class{1,1}==1);
sub_abnormal=find(X_rem.class{1,1}==2);

%% univariate analysis
s=size(T0);
for k=1:s(2)
    [h, p_meta(k), ~, tt]=ttest2(T0(sub_normal,k), T0(sub_abnormal,k), 'alpha', 0.05, 'vartype','unequal');
end

%% preprocessing
% centering
T0_centered=T0-repmat(mean(T0,1),size(T0,1),1);
T0_T0=T0_centered;
% scale
for i=1:size(T0_T0,2)
    rms=sqrt(mean(T0_T0(:,i).^2));
    T0_scal(:,i)=T0_T0(:,i)/rms;
end


%% perform svd
[U,S,V] = svd(T0_scal);

%% check about the p-val in each component for the subjects mode
nm_comp=6;
for r=1:nm_comp
    r1=U(sub_normal,r);
    r2=U(sub_abnormal,r);
    [h, p(r), ~, tt]=ttest2(r1, r2, 'alpha', 0.05, 'vartype','unequal');
end
pmin=min(p)
index_pmin=find(p==pmin)

%% plot the explained variance
Sdiag=diag(S);
Epval=cumsum(Sdiag(1:length(Sdiag)).^2)/sum(Sdiag.^2)*100;
figure
plot(1:length(Epval), Epval,'--*','Markersize',8,'LineWidth',2)
grid on
xlabel('number of components')
ylabel('explained variance')
title('PCA model')
set(gca,'Fontsize',15)

%%
PC_exp(1)=Epval(1);
PC_exp(2:length(Epval))=Epval(2:end)-Epval(1:end-1);
PC_exp=round(PC_exp,0);
%%



% %% plot the factor plots
labelss = {'Ins','GLC','Pyr','Lac','Ala','Bhb'};
k=0;
% x_points=1:size(U,1);
% figure
% for i=1:nm_comp
%     subplot(2,nm_comp,i)
%     plot(x_points(sub_normal), U(sub_normal,i),'rs','Markersize',6,'LineWidth',2);
%     
%     hold on
%     plot(x_points(sub_abnormal), U(sub_abnormal,i),'bs','Markersize',6,'LineWidth',2);
%     grid on
%     ylabel(['PC',num2str(i),'-',num2str(PC_exp(i)),'%'])
%     set(gca,'Fontsize',13)
% end
% for i=1:nm_comp
%     subplot(2,nm_comp,nm_comp+i)
%     plot(1:length(V(:,i)), V(:,i),'-o','LineWidth',2);
%     grid on
%     ylabel(['PC',num2str(i),'-',num2str(PC_exp(i)),'%'])
%     xtickangle(90)
%     set(gca,'xtick',1:length(V(:,i)),'xticklabel',labelss)
%     set(gca,'Fontsize',13)
% end




%% the sactter plot
labelss = {'Ins','GLC','Pyr','Lac','Ala','Bhb'};
k=0;
for i=1:nm_comp
    for j=i+1:nm_comp
        figure
        subplot(2,1,1)
        plot(U(sub_normal,i), U(sub_normal,j),'ro','MarkerFaceColor','r');
        hold on
        plot(U(sub_abnormal,i), U(sub_abnormal,j),'bo','MarkerFaceColor','b');
        
        xlabel(['PC',num2str(i),'-',num2str(PC_exp(i)),'%'])
        ylabel(['PC',num2str(j),'-',num2str(PC_exp(j)),'%'])
        
        set(gca,'Fontsize',13)
        subplot(2,1,2)
        plot(V(:,i), V(:,j),'ro','MarkerFaceColor','r');
        text(V(:,i), V(:,j),labelss,'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize', 15,'fontweight', 'bold')
        xlabel(['PC',num2str(i),'-',num2str(PC_exp(i)),'%'])
        ylabel(['PC',num2str(j),'-',num2str(PC_exp(j)),'%'])
        grid on
        set(gca,'Fontsize',13)
    end
end

