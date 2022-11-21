% This script shows how to fit the PCA model using the svd algorithm to the
% fasting-state (T0) real data (only male subjects and 162 measurements) with the missing
% values replaced by the approximations



clear all
clc
close all
%%
load('NMR+Gene+SNP+Metav3_July4_2022.mat') %% load data

%% remove outliers
pid_list = str2num(NMR.label{1});
outlier_pid =[142,79,342,250,90,335,312];
NMR_remove=remove_outliers(NMR,outlier_pid);

%% remove the metabolites with lots of missing values
NMR_remove=removeisnan(NMR_remove);

%% select the large set of metabolites (162 measurements)
selected_metabolite_index = find(NMR_remove.class{3,1}==1 | NMR_remove.class{3,1}==2);
NMR_remove= NMR_remove(:,:,selected_metabolite_index);

%% choose male(2)/female(1)
selected_gender_index =find(NMR_remove.class{1,1}==2);
NMR_remove= NMR_remove(selected_gender_index,:,:);

%% get the raw fasting-state data
X_slect = squeeze(NMR_remove.data(:,1,:));

%% replace the missing values with the approximated values by spca
load('data_T0_approx.mat','data_T0_approx')
X_slect_new=X_slect;
[row,col] = find(isnan(X_slect));
for i=1:length(row)
    X_slect_new(row(i),col(i))=data_T0_approx(row(i),col(i));
end
%% preprocessing
s=size(X_slect_new);
X_center=X_slect_new-repmat(mean(X_slect_new,1),s(1),1);
for i=1:s(2)
    temp=X_center(:,i);
    X_pre(:,i)=temp/sqrt(mean(temp.^2));
end
X=X_pre;

%%  PCA using svd
nm_comp=2;
[U, S, V]=svd(X);
Fac=ktensor(diag(S(1:nm_comp,1:nm_comp)),U(:,1:nm_comp), V(:,1:nm_comp));
Fac=normalize(Fac);

%% compute the explained variance
sigma=diag(S);
expl=cumsum(sigma.^2)/sum(sigma.^2);
expl=100*expl;
Expl_comp(1)=expl(1);
Expl_comp(2:length(expl))=expl(2:end)-expl(1:end-1);
Expl_comp=round(Expl_comp,0);

%% Flip the factors for visulization
Fac.U{1}(:,1)=-Fac.U{1}(:,1);Fac.U{2}(:,1)=-Fac.U{2}(:,1);
Fac.U{1}(:,2)=-Fac.U{1}(:,2);Fac.U{2}(:,2)=-Fac.U{2}(:,2);

%% BMI index
index_Under=find(NMR_remove.class{1,2}==4);
index_normal=find(NMR_remove.class{1,2}==1);
index_obesity=find(NMR_remove.class{1,2}==2);
index_over=find(NMR_remove.class{1,2}==3);
index_Under_normal=[index_Under,index_normal];
index_over_obesity=[index_over,index_obesity];

%% BMI separation and boxplot of the separation component
for i=1:nm_comp
    [~,p_val(i)] = ttest2(Fac.U{1}(index_Under_normal,i),Fac.U{1}(index_over_obesity,i),'Vartype','unequal');
end

x1=Fac.U{1}(index_Under_normal,1);
x2=Fac.U{1}(index_over_obesity,1);
x=[x1; x2];
g1 = repmat({'High BMI'},length(x1),1);
g2 = repmat({'Low BMI'},length(x2),1);
g = [g1; g2];
boxplot(x,g)
ylabel('PC1')
set(gca,'Fontsize',16)


%% metabolites index
plot_color = ['g','c','m','k'];
plot_shape = ['s','s','s','s'];
index_HDL=find(NMR_remove.class{3,3}==1);
index_LDL=find(NMR_remove.class{3,3}==2);
index_VLDL=find(NMR_remove.class{3,3}==4);
index_Rest=setdiff(1:length(Fac.U{2}(:,1)), [index_HDL, index_LDL, index_VLDL]);

%% scatter plot by PCs
Mark_value=5;
subplot(1,2,1)
for i=1:nm_comp
    for j=i+1:nm_comp
        p(1)=plot(Fac.U{1}(index_Under_normal,i),Fac.U{1}(index_Under_normal,j),'rs','MarkerSize',Mark_value,...
            'MarkerEdgeColor','r','MarkerFaceColor','r');
        hold on
        p(2)=plot(Fac.U{1}(index_over_obesity,i),Fac.U{1}(index_over_obesity,j),'bs','MarkerSize',Mark_value,...
            'MarkerEdgeColor','b','MarkerFaceColor','b');
        hold on
        xlim_value=[min([min(Fac.U{1}(:,i)), min(Fac.U{1}(:,j))]),...
            max([max(Fac.U{1}(:,i)), max(Fac.U{1}(:,j))])];
        xlim(xlim_value)
        ylim(xlim_value)
        p(3)=plot(xlim_value, [0 0],'k-');
        hold on
        p(4)=plot([0 0], xlim_value,'k-');
        legend(p(1:2), 'Low BMI','High BMI')
        xlabel(['PC' num2str(i),'-',num2str(Expl_comp(i)),'%'])
        ylabel(['PC' num2str(j),'-',num2str(Expl_comp(j)),'%'])
        title('subjects')
        set(gca,'Fontsize',16)
    end
end
clear p  xlim_value
subplot(1,2,2)
for i=1:nm_comp
    for j = i+1:nm_comp
        x_value =Fac.U{2}(:,i);
        y_value =Fac.U{2}(:,j);
        ii=1;p(ii)=plot(x_value(index_HDL),y_value(index_HDL),[plot_color(ii) plot_shape(ii)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(ii)]);
        hold on
        ii=2;p(ii)=plot(x_value(index_LDL),y_value(index_LDL),[plot_color(ii) plot_shape(ii)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(ii)]);
        hold on
        ii=3;p(ii)=plot(x_value(index_VLDL),y_value(index_VLDL),[plot_color(ii) plot_shape(ii)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(ii)]);
        hold on
        ii=4;p(ii)=plot(x_value(index_Rest),y_value(index_Rest),[plot_color(ii) plot_shape(ii)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(ii)]);
        hold on
        xlim_value=[min([min(Fac.U{2}(:,i)), min(Fac.U{2}(:,j))]),...
            max([max(Fac.U{2}(:,i)), max(Fac.U{2}(:,j))])];
        xlim(xlim_value)
        ylim(xlim_value)
        ii=5;p(ii)=plot(xlim_value, [0 0],'k-');
        hold on
        ii=6;p(ii)=plot([0 0], xlim_value,'k-');
        hold on
        legend(p(1:4), 'HDL','LDL','VLDL','Rest')
        xlabel(['PC' num2str(i),'-',num2str(Expl_comp(i)),'%'])
        ylabel(['PC' num2str(j),'-',num2str(Expl_comp(j)),'%'])
        title('metabolites')
        xticks ([-0.08 0  0.08, 0.15])
        yticks ([-0.08 0  0.08, 0.15])
        set(gca,'Fontsize',16)
    end
    
end



