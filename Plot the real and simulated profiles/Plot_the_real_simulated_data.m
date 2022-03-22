%%
clear all
clc
close all
% addpath('/Users/ll/Documents/MATLAB/toolbox/PLS_Toolbox_891') % for dataset


%% load the data
load('NMR+Gene+SNP+Metav3_Feb14_2022.mat')
%% pick out the common metabolites with the simulated data
% The intersected metabolites are {'Ins','GLC','Pyr','Lac','Ala','Bhb'};
a=NMR.label{3,1};
b=['Insulin        ';'Glucose        ';'Pyruvate       ';'Lactate        ';...
    'Ala            ';'Glycerol       ';'bOHbutyrate    ';
    'Total-TG       ';'Total-C        ';   'Unsaturation   '];
% b=['Insulin        ';'Glucose        ';'Pyruvate       ';'Lactate        ';...
%     'Ala            ';'Glycerol       ';'bOHbutyrate    ';
%     'Total-TG       ';'Total-C        ';'Total-FA       '];

% do not include Glycerol since most of them are missing. Missing check
% will remove this metabolite
selected_index=find( ismember(a, b, 'rows'));
NMR=NMR(:,:,selected_index); % If you use other metabolites,just remove this line

%% remove subjects and metabolites that have too many nans
NMR=removeisnan(NMR);  

%% reorder the metabolites as {'Ins','GLC','Pyr','Lac','Ala','Bhb'};
index=[9 5 7 6 4 8 3 2 1 ];
NMR=NMR(:,:,index);
labelss = {'Ins','GLC','Pyr','Lac','Ala','Bhb','FFA','TG','Chol'};

%% parameters for selecting the data
intersec_gene=1;%% 0 means not intersected data;1 means intersected data
gender='male'; %% options are 'female','male','full'
meta_type=0; %% options are 1,2 and others;1 means Age_selected(small) sets;2 means Gozde select(large) sets.
data_type='full dynamic'; %% options are: 'T0', 'T0 corrected', 'full dynamic'

%% select data (gene-intersection,male/female/full) from the real data
NMR=select_data(NMR,Gene,intersec_gene,gender,meta_type);

%% consider the full dynamic data
X_slect=NMR.data;

%% BMI index
index_Under=find(NMR.class{1,2}==4);
index_normal=find(NMR.class{1,2}==1);
index_obesity=find(NMR.class{1,2}==2);
index_over=find(NMR.class{1,2}==3);
index_Under_normal=[index_Under,index_normal];
index_over_obesity=[index_over,index_obesity];

%% load the simulated data and pick out the common metabolites 
load('Simu_10meta_8time4hour_beta03.mat')
index=[1 2 3 4 5 7 8 9 10];
X_orig=X_orig(:,index,:);
X_orig.data(:,1,:)=X_orig.data(:,1,:)*10^9; % Insulin in the real data is pM while it is mM in the simulated data
%% find index for normal and abnormal subjects in simulated data
 sub_normal=find(X_orig.class{1,1}==2);
sub_abnormal=find(X_orig.class{1,1}==1);



%% plot raw data
% timelabel=str2num(NMR.label{2,1});
figure
k=0;
for i=1:size(X_slect,3)
    k=k+1;
    subplot(3,4,k)
    plot(str2num(NMR.label{2,1}),squeeze(nanmean(X_slect(index_normal,:,i),1)),'rs','Linewidth',3)
    hold on
    plot(str2num(NMR.label{2,1}),squeeze(nanmean(X_orig.data( sub_normal,i,:),1)),'-r','Linewidth',3)
    xticks([0:1:4])
    xlabel('hour')
    ylabel('concentration')
    title([labelss(i)])
    grid on
    set(gca,'Fontsize',14)
    real_normal_start(i)=squeeze(nanmean(X_slect(index_normal,1,i),1));
end
 legend('real','simulated')



for i=1:9
    figure
    subplot(1,2,1)
    plot(str2num(NMR.label{2,1}),squeeze(nanmean(X_slect(index_normal,:,i),1)),'--rs','Linewidth',3)
     xticks([0:1:4])
    title([labelss(i),'-real'])
    set(gca,'Fontsize',14)
     xlabel('hour')
    ylabel('concentration')
    subplot(1,2,2)
    plot(str2num(NMR.label{2,1}),squeeze(nanmean(X_orig.data( sub_normal,i,:),1)),'-rs','Linewidth',3)
     xticks([0:1:4])
    title([labelss(i),'-simulated'])
    xlabel('hour')
    ylabel('concentration')
    set(gca,'Fontsize',14)
end







