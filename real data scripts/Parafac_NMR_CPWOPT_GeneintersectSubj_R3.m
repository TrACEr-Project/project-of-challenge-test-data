clear all
clc
close all
% need to add tensor toolbox 
% addpath('/Users/ll/Documents/MATLAB/toolbox/tensor_toolbox_v3.1')
% addpath('/Users/ll/Documents/MATLAB/toolbox/PLS_Toolbox_891') % for use of datset


%% load the real data
load('NMR+Gene+SNP+Metav3_Feb14_2022.mat')

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

%% consider the full data
X_slect=NMR.data;

%% BMI index
index_Under=find(NMR.class{1,2}==4);
index_normal=find(NMR.class{1,2}==1);
index_obesity=find(NMR.class{1,2}==2);
index_over=find(NMR.class{1,2}==3);
index_Under_normal=[index_Under,index_normal];
index_over_obesity=[index_over,index_obesity];

%% plot raw data

k=0;
figure
for i=1:6
    k=k+1;
    subplot(2,3,k)
    for j=index_Under_normal
        plot(squeeze(X_slect(j,:,i)),'r')
        hold on
    end
    for j=index_over_obesity
        plot(squeeze(X_slect(j,:,i)),'b')
        hold on
    end
    plot(squeeze(nanmean(X_slect(index_Under_normal,:,i),1)),'k','Linewidth',3)
    hold on
    plot(squeeze(nanmean(X_slect(index_over_obesity,:,i),1)),'g','Linewidth',3)
    xlim([1, 8])
     title(labelss(i))
end

%% preprocessing
X=tensor(preprocess(X_slect));

%% plot preprocessed data
k=0;
figure
for i=1:6
    k=k+1;
    subplot(2,3,k)
    for j=index_Under_normal
        plot(squeeze(X.data(j,:,i)),'r')
        hold on
    end
    for j=index_over_obesity
        plot(squeeze(X.data(j,:,i)),'b')
        hold on
    end
    plot(squeeze(nanmean(X.data(index_Under_normal,:,i),1)),'k','Linewidth',3)
    hold on
    plot(squeeze(nanmean(X.data(index_over_obesity,:,i),1)),'g','Linewidth',3)
    xlim([1, 8])
     title(labelss(i))
end



%%  cp_wopt model
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
%% uniqueness test
% 0 -> NOT unique
% 1 -> Unique
% 2 -> Inconclusive, need more random starts
good_flag = find(goodness_X1(:) == 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL.' | goodness_X1(:) == 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH.');
% good_flag =1:nb_starts;
if length(good_flag)>=1
    F_round = round(goodness_X(good_flag,2),6);
    best_F_index = good_flag(F_round == min(F_round));
    if length(best_F_index) < 2 %Try for 1e-8, but can do 1e-7
        F_round = round(goodness_X(good_flag,2),0); %Round F to 1e-7,TOO TIGHT?
        best_F_index = good_flag(F_round == min(F_round));%Finds best F values
    end
    %     flag = 0;
else
    F_round = round(goodness_X(:,2),10);
    best_F_index = find(F_round == min(F_round));
    %     flag = 1;
end

eps = .05; %Arbtitraly picked, ideas for a values are appreciated
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
Fac = Fac_X{index(1)};
uniqueness=unique_test
fit=   goodness_X(index(1),1)
exit_info=goodness_X1(best_F_index(1))
Consistency = corcond(X.data,normalize(Fac,1),[],0)
tc=TC(Fac.U)







%% Explore BMI separation

% compute the p-val
for r=1:nm_comp
    r1=Fac.U{1}(index_Under_normal,r);
    r2=Fac.U{1}(index_over_obesity,r);
    [h, p(r), ~, tt]=ttest2(r1, r2, 'alpha', 0.05, 'vartype','unequal');
end
[pmin,index_pmin]=min(p);
% plot the subjects mode 
xvalue=1:length(Fac.U{1}(:,1));
figure
for i=1:nm_comp
    subplot(nm_comp,1,i)
    plot(xvalue(index_Under_normal),Fac.U{1}(index_Under_normal,i),'rs','MarkerSize',8,...
        'MarkerEdgeColor','r','MarkerFaceColor','r')
    hold on
    plot(xvalue(index_over_obesity),Fac.U{1}(index_over_obesity,i),'bs','MarkerSize',8,...
        'MarkerEdgeColor','b','MarkerFaceColor','b')
    set(gca,'FontSize', 18)
    legend('Under+normal','over+obesity')
     ylabel(['Comp',num2str(i)])
end

% scatter plot of the subjects mode
for i=1:nm_comp
    for j=i+1:nm_comp
    figure
    plot(Fac.U{1}(index_Under_normal,i),Fac.U{1}(index_Under_normal,j),'rs','MarkerSize',8,...
        'MarkerEdgeColor','r','MarkerFaceColor','r')
    hold on
    plot(Fac.U{1}(index_over_obesity,i),Fac.U{1}(index_over_obesity,j),'bs','MarkerSize',8,...
        'MarkerEdgeColor','b','MarkerFaceColor','b')
    set(gca,'FontSize', 18)
    legend('Under+normal','over+obesity')
    xlabel(['Comp',num2str(i)])
    ylabel(['Comp',num2str(j)])
    end
end


%% plot the factors
Z2={'subjects','time','metabolites'};
figure
nm_comp=3;
for i=1:3
    if i==3
    for j=1:nm_comp
    subplot(3,nm_comp,(i-1)*nm_comp+j)
    plot(Fac.U{i}(:,j),'-s','Markersize',6,'LineWidth',2)
    set(gca,'xtick',1:length(Fac.U{i}(:,j)),'xticklabel',labelss)
    xtickangle(90)
    ylabel(['comp',num2str(j)])
    ylim([min(min(Fac.U{i})), max(max(Fac.U{i}))])
    grid on
    end
    elseif i==1
         for j=1:nm_comp
      subplot(3,nm_comp,(i-1)*nm_comp+j)
      plot(xvalue(index_Under_normal),Fac.U{i}(index_Under_normal,j),'rs','MarkerSize',6,...
        'MarkerEdgeColor','r','MarkerFaceColor','r')
    hold on
    plot(xvalue(index_over_obesity),Fac.U{i}(index_over_obesity,j),'bs','MarkerSize',6,...
        'MarkerEdgeColor','b','MarkerFaceColor','b')
    ylabel(['comp',num2str(j)])
    ylim([min(min(Fac.U{i})), max(max(Fac.U{i}))])
    grid on
         end  
    else
    for j=1:nm_comp
    subplot(3,nm_comp,(i-1)*nm_comp+j)
    plot(1:length(Fac.U{i}(:,j)),Fac.U{i}(:,j),'-s','Markersize',6,'LineWidth',2)
    
        xlabel(Z2{i})
    
   ylabel(['comp',num2str(j)])
   ylim([min(min(Fac.U{i})), max(max(Fac.U{i}))])
    grid on
    end
    end
end

    
% %%
% check_outliers_CP(X,Fac)
