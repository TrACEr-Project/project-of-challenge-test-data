% This script shows an example to fit a CP model to the full-dynamic data from the simulated datasets.

% We use Tensor Toolbox as well as the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C
% In addition, parts of the scripts may require the dataset object (https://eigenvector.com/software/dataset-object/), publically available.
% For core consistency computation, we also use the corcord function from the Nway toolbox (http://www.models.life.ku.dk/nwaytoolbox).


%% add auxilary functions to path
addpath(genpath('.\functions'))
%% add dataset path
addpath(genpath('..\simulated_datasets\Betacell_dysfunction'))
%% add other apckages to your path!
addpath(genpath('...\tensor_toolbox-v3.1')) %Tensor toolbox is needed;  MATLAB Tensor Toolbox. Copyright 2017, Sandia Corporation, http://www.tensortoolbox.org/
addpath(genpath('...\L-BFGS-B-C-master')) % LBFGS-B implementation is needed; download here: https://github.com/stephenbeckr/L-BFGS-B-C
addpath(genpath('...\nway331')) % Nway toolbox is needed for computing core consistency; download here: http://www.models.life.ku.dk/nwaytoolbox
addpath(genpath('...\dataset')) % dataset object is needed; download here: https://eigenvector.com/software/dataset-object/

%%  load dataset
load('Simu_6meta_8time_alpha02_betacell_balance.mat','X_orig')


%% remove subjects with blow-up solution when solving ODE / remove outliers
nr_sub_zeros=find(X_orig.class{1,2}==2); % subjects with blow-up solution
pid_list=str2num(X_orig.label{1});
outlier_index=[nr_sub_zeros];
outlier_pid=pid_list(outlier_index);
X_rem=removesubject(X_orig,outlier_pid);


%% full dynamic data
X=X_rem.data;
s=size(X);


%% find index for normal and abnormal subjects
sub_normal=find(X_rem.class{1,1}==1);
sub_abnormal=find(X_rem.class{1,1}==2);


%% plot the raw data
labelss = {'Ins','GLC','Pyr','Lac','Ala','Bhb'};
time_value=[0 0.25 0.5 1 1.5 2 2.5 4];
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
unique_test = unique_test_CP(Fac_X, goodness_X(:,2), goodness_X1)

%%  Get the best Factorization & fit,TC, CC info
[er,index]=sort(goodness_X(:,2),'ascend');
Fac = Fac_X{index(1)};
out_X_best = out_X{index(1)};
exit_info=goodness_X1(index(1))
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
time_value=[0 0.25 0.5 1 1.5 2 2.5 4];
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
            ylabel(['b_',num2str(j)],'fontweight','bold')
            set(gca,'Fontsize',15)
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
            if j==4
                legend('control','diseased')
            end
            ylim([min(min(Fac.U{i})),max(max(Fac.U{i}))])
            ylabel(['a_',num2str(j)],'fontweight','bold')
            xlabel(Z2{i})
            set(gca,'Fontsize',15)
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
            ylabel(['c_',num2str(j)],'fontweight','bold')
            xlabel('time(h)')
            set(gca,'Fontsize',15)
        end
    end
    
    
end



