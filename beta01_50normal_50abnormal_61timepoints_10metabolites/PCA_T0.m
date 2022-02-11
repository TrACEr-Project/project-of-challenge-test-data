clear all
clc
% close all

%% load data

load ('Simu_10meta_61time_beta01.mat','X_orig');
labelss = {'Ins_B','GLC_B','Pyr_B','Lac_B','Ala_B','Glyc_B','Bhb_B','FFA_B',...
    'TG_B','Chol_B'};
%% remove outliers
pid_list=str2num(X_orig.label{1});
outlier_index=[];
outlier_pid=pid_list(outlier_index);
X_rem=removeoutlier(X_orig,outlier_pid);
X=X_rem.data;

%% find index for normal and abnormal subjects
sub_normal=find(X_rem.class{1,1}==1);
sub_abnormal=find(X_rem.class{1,1}==2);


%% preprocessing
T0=X(:,:,1);
% centering
T0_centered=T0-repmat(mean(T0,1),size(T0,1),1);
T0_T0=T0_centered;
% scale
for i=1:size(T0_T0,2)
    rms=sqrt(mean(T0_T0(:,i).^2));
    T0_scal(:,i)=T0_T0(:,i)/rms;
end
s=size(T0_scal);
%% plot the raw data
figure
for i=1:s(2)
    subplot(4,3,i)
    for j=sub_normal
        plot(squeeze(T0(j,i)),'ro')
        hold on
    end
    for j=sub_abnormal
        plot(squeeze(T0(j,i)),'bo')
        hold on
    end
    title(labelss(i))
    set(gca,'FontSize', 15)
end

%% plot the preprocessed data
figure
for i=1:s(2)
    subplot(4,3,i)
    for j=sub_normal
        plot(squeeze(T0_scal(j,i)),'ro')
        hold on
    end
    for j=sub_abnormal
        plot(squeeze(T0_scal(j,i)),'bo')
        hold on
    end
    title(labelss(i))
    set(gca,'FontSize', 15)
end

%%
[U,S,V] = svd(T0_scal);


%%
Fac.lambda=diag(S);
Fac.U{1}=U;
Fac.U{2}=V;
Fac_PCA_T0=Fac;

nm_comp=10;
for r=1:nm_comp
    r1=U(sub_normal,r);
    r2=U(sub_abnormal,r);
    [h, p(r), ~, tt]=ttest2(r1, r2, 'alpha', 0.05, 'vartype','unequal');
end
[pmin,index_pmin]=min(p)
%%
figure
nm_comp1=5;
for i=1:nm_comp1
    subplot(2,nm_comp1,i)
    plot(Fac_PCA_T0.U{1}(:,i),'-o','LineWidth',2);
    hold on
    ylabel(['comp', num2str(i) ],'FontSize', 10)
end
for i=1:nm_comp1
    subplot(2,nm_comp1,nm_comp1+i)
    plot(Fac_PCA_T0.U{2}(:,i),'-o','LineWidth',2.4)
    hold on
    ylabel(['comp', num2str(i) ],'FontSize', 10)
    grid on
    set(gca,'xtick',1:length(Fac.U{2}(:,i)),'xticklabel',labelss)
    xtickangle(90)
    
end

%%
figure
nm_comp2=5;
for i=1:nm_comp2
    subplot(2,nm_comp2,i)
    plot(Fac_PCA_T0.U{1}(:,nm_comp1+i),'-o','LineWidth',2);
    hold on
    ylabel(['comp', num2str(nm_comp1+i) ],'FontSize', 10)
end
for i=1:nm_comp1
    subplot(2,nm_comp2,nm_comp2+i)
    plot(Fac_PCA_T0.U{2}(:,nm_comp1+i),'-o','LineWidth',2.4)
    hold on
    ylabel(['comp', num2str(nm_comp1+i) ],'FontSize', 10)
    grid on
    set(gca,'xtick',1:length(Fac.U{2}(:,nm_comp1+i)),'xticklabel',labelss)
    xtickangle(90)
    
end



%%

Sdiag=diag(S);
Epval=cumsum(Sdiag(1:length(Sdiag)).^2)/sum(Sdiag.^2)*100;
figure
plot(1:length(Epval), Epval,'--*','Markersize',8,'LineWidth',2)
xlabel('number of components')
ylabel('explained variance')
title('PCA model')
set(gca,'Fontsize',15)

