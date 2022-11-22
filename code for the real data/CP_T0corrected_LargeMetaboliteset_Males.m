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
load('NMR+Gene+SNP+Metav3_July4_2022.mat') %% load the real data

%% remove outliers
pid_list = str2num(NMR.label{1});
outlier_pid =[142,79,342,250,90,335,312];
NMR_remove=removesubject(NMR,outlier_pid);

%% remove the metabolites with lots of missing values
NMR_remove=removeisnan(NMR_remove);

%% select the large set of metabolites (162 measurements)
selected_metabolite_index = find(NMR_remove.class{3,1}==1 | NMR_remove.class{3,1}==2);
NMR_remove= NMR_remove(:,:,selected_metabolite_index);

%% choose male(2)/female(1)
selected_gender_index =find(NMR_remove.class{1,1}==2);
NMR_remove= NMR_remove(selected_gender_index,:,:);


%% Get the raw T0 corrected data
X_slect = NMR_remove.data(:,2:end,:)-NMR_remove.data(:,1,:);

%% preprocessing
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
nb_starts =60;
nm_comp=2;
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

%% BMI index
index_Under=find(NMR_remove.class{1,2}==4);
index_normal=find(NMR_remove.class{1,2}==1);
index_obesity=find(NMR_remove.class{1,2}==2);
index_over=find(NMR_remove.class{1,2}==3);
index_Under_normal=[index_Under,index_normal];
index_over_obesity=[index_over,index_obesity];

%% boxplot of the BMI separation component
for i=1:nm_comp
    [~,p_val(i)] = ttest2(Fac.U{1}(index_Under_normal,i),Fac.U{1}(index_over_obesity,i),'Vartype','unequal');
end

x1=Fac.U{1}(index_Under_normal,2);
x2=Fac.U{1}(index_over_obesity,2);
x=[x1; x2];
g1 = repmat({'High BMI'},length(x1),1);
g2 = repmat({'Low BMI'},length(x2),1);
g = [g1; g2];
boxplot(x,g)
ylabel(['a_',num2str(2)],'fontweight','bold')
set(gca,'Fontsize',16)


%% metabolites index
plot_color = ['g','c','m','k'];
plot_shape = ['s','s','s','s'];
index_HDL=find(NMR_remove.class{3,3}==1);
index_LDL=find(NMR_remove.class{3,3}==2);
index_VLDL=find(NMR_remove.class{3,3}==4);
index_Rest=setdiff(1:length(Fac.U{2}(:,1)), [index_HDL, index_LDL, index_VLDL]);

%% plot of the CP factors
f=figure;
Mark_value=5;
time_values=[0.25 0.5 1 1.5 2 2.5 4];
k=0;
for ii=1:3
    for jj = 1:nm_comp
        k=k+1;
        subplot(3, nm_comp, k)
        x_value = 1:length(Fac.U{ii}(:,1));
        y_value =Fac.U{ii}(:,jj);
        if ii==2
            
            i=1;p(1)=plot(x_value(index_HDL),y_value(index_HDL),[plot_color(i) plot_shape(i)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(i)]);
            hold on
            i=2;p(2)=plot(x_value(index_LDL),y_value(index_LDL),[plot_color(i) plot_shape(i)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(i)]);
            hold on
            i=3;p(3)=plot(x_value(index_VLDL),y_value(index_VLDL),[plot_color(i) plot_shape(i)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(i)]);
            hold on
            i=4;p(4)=plot(x_value(index_Rest),y_value(index_Rest),[plot_color(i) plot_shape(i)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(i)]);
            hold on
            p(5)=plot([min(x_value) max(x_value)],[0.08 0.08],'--k','Linewidth',0.8);
            hold on
            p(6)=plot([min(x_value) max(x_value)],[-0.08 -0.08],'--k','Linewidth',0.8);
            hold on
            p(7)=plot([min(x_value) max(x_value)],[0 0],'--k','Linewidth',0.8);
            xlabel(['metabolites'])
            ylabel(['b_',num2str(jj)],'fontweight','bold')
            legend([p(1) p(2) p(3) p(4)],'HDL','LDL','VLDL','Rest','Fontsize',12)
            ylim([min(min(Fac.U{ii})),max(max(Fac.U{ii}))])
            xlim([min(x_value) max(x_value)])
            yticks([-0.15 -0.08 0 0.08 0.15])
            set(gca,'Fontsize',18)
        elseif ii==1
            
            plot(x_value(index_Under_normal),Fac.U{ii}(index_Under_normal,jj),'rs','MarkerSize',Mark_value,...
                'MarkerEdgeColor','r','MarkerFaceColor','r')
            hold on
            plot(x_value(index_over_obesity),Fac.U{ii}(index_over_obesity,jj),'bs','MarkerSize',Mark_value,...
                'MarkerEdgeColor','b','MarkerFaceColor','b')
            set(gca,'Fontsize',18)
            if jj==2
                legend('Low BMI','High BMI','Location','northeast','Fontsize',12)
            end
            ylabel(['a_',num2str(jj)],'fontweight','bold')
            xlabel(['subjects'])
            ylim([min(min(Fac.U{ii})),max(max(Fac.U{ii}))])
            xlim([min(x_value) max(x_value)])
            set(gca,'Fontsize',18)
        else
            plot(time_values,Fac.U{ii}(:,jj),'-*','Linewidth',3)
            xticks(0:4)
            ylabel(['c_',num2str(jj)],'fontweight','bold')
            xlabel(['time(h)'])
            ylim([min(min(Fac.U{ii})),max(max(Fac.U{ii}))])
            set(gca,'Fontsize',18)
        end
        
    end
end




%% scatter plot of the subjects and metabolites
Mark_value=8;
f=figure;
subplot(1,2,1)
for i=1:nm_comp
    for j=i+1:nm_comp
        p(1)=plot(Fac.U{1}(index_Under_normal,i),Fac.U{1}(index_Under_normal,j),'rs','MarkerSize',Mark_value,...
            'MarkerEdgeColor','r','MarkerFaceColor','r');
        hold on
        p(2)=plot(Fac.U{1}(index_over_obesity,i),Fac.U{1}(index_over_obesity,j),'bs','MarkerSize',Mark_value,...
            'MarkerEdgeColor','b','MarkerFaceColor','b');
        hold on
        plot([min(min(Fac.U{1})),max(max(Fac.U{1})) ], [0 0]);
        hold on
        plot( [0 0], [min(min(Fac.U{1})),max(max(Fac.U{1})) ]);
        legend([p(1) p(2)],'Low BMI','High BMI','Fontsize',14)
        xlim([min(min(Fac.U{1})),max(max(Fac.U{1})) ])
        ylim([min(min(Fac.U{1})),max(max(Fac.U{1})) ])
        xlabel(['a_',num2str(i)],'fontweight','bold')
        ylabel(['a_',num2str(j)],'fontweight','bold')
        title('subjects')
        set(gca,'FontSize', 16)
    end
end
subplot(1,2,2)
for cpi=1:nm_comp
    for cpj = cpi+1:nm_comp
        x_value = Fac.U{2}(:,cpi);
        y_value =Fac.U{2}(:,cpj);
        i=1;p(i)=plot(x_value(index_HDL),y_value(index_HDL),[plot_color(i) plot_shape(i)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(i)]);
        hold on
        i=2;p(i)=plot(x_value(index_LDL),y_value(index_LDL),[plot_color(i) plot_shape(i)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(i)]);
        hold on
        i=3;p(i)=plot(x_value(index_Rest),y_value(index_Rest),[plot_color(i) plot_shape(i)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(i)]);
        hold on
        i=4;p(i)=plot(x_value(index_VLDL),y_value(index_VLDL),[plot_color(i) plot_shape(i)],'MarkerSize',Mark_value,'MarkerFaceColor',[plot_color(i)]);
        hold on
        plot([min(min(Fac.U{2})),max(max(Fac.U{2})) ], [0 0]);
        hold on
        plot( [0 0], [min(min(Fac.U{2})),max(max(Fac.U{2})) ]);
        xlim([min(min(Fac.U{2})),max(max(Fac.U{2})) ])
        ylim([min(min(Fac.U{2})),max(max(Fac.U{2})) ])
        xlabel(['b_',num2str(cpi)],'fontweight','bold')
        ylabel(['b_',num2str(cpj)],'fontweight','bold')
        legend([p(1) p(2) p(3) p(4)],'HDL','LDL','VLDL','Rest','Fontsize',14)
        title('metabolites')
        set(gca,'Fontsize',16)
    end
    
end
