%%This is for selecting data we want to use from the real data

function [NMR_select]=select_data(NMR,Gene,intersec_gene,gender,meta_type)


%% use gene and NMR intersection or not
switch intersec_gene
    case 0
        NMR_select=NMR;
    case 1
        Gene_X = Gene;
        Gene_X = Gene_X(:,[1:12 16:end]);%% gene data needs to remove a threshold.
        pid_NMR = str2num(NMR.label{1});
        pid_Gene = str2num(Gene_X.label{1,1});
        [common,subid_NMR,subid_Gene] = intersect(pid_NMR,pid_Gene);
        NMR_select=NMR(subid_NMR,:,:);
       

end

%%  select metabolites type: age_select 1 , gozde select 2
switch meta_type
    case 2
        selected_index = find(NMR_select.class{3,4}>0);
    case 1
        selected_index = find(NMR_select.class{3,1}>0);
    otherwise
        selected_index = 1:size(NMR_select.data,3);
        disp(['all metabolites']);
end
NMR_select=NMR_select(:,:,selected_index);
%% choose male/female/full
switch gender
    case  'full gender'
        selected_index = 1:size(NMR_select.data,1);
    case 'female'
        selected_index = find(NMR_select.class{1,1}==1);
        
    case  'male'
        selected_index =find(NMR_select.class{1,1}==2);
        
end

NMR_select=NMR_select(selected_index,:,:);
% Metainfo=Metainfo(selected_index,:);

end

%% function for remove subjects/metabolites with lots of nans
function X= removeisnan(X)
% check if any subject has no data, or any metabolite has no data.
% need a loop
do_loop = 1;

while do_loop
    % check subject
    subject_id = 1:size(X,1);
    kept_subject = ones(size(X,1),1);
    do_loop = 0;
    for pid = 1:size(X,1)
        
        data_count = sum(~isnan(X.data(pid,:,:)),'all');
        if data_count<0.3*size(X,2)*size(X,3) % 0 means less data in this pid
            disp(['Pid ' num2str(pid) 'contains less data.'])
            kept_subject(pid) = 0;
            do_loop = 1;
        end
    end
    kept_subject_id = subject_id(find(kept_subject==1));
    X = X(kept_subject_id,:,:);
    % check metabolite
    metabolite_id = 1:size(X,3);
    kept_metabolite = ones(size(X,3),1);
    
    for metaid = 1:size(X,3)
        
        data_count = sum(~isnan(X.data(:,:,metaid)),'all');
        %disp(['Meta Num ' num2str(metaid) ' ' X.label{3,1}(metaid,:) ' has ' num2str(data_count) ' data in ' ])
        if data_count<0.3*size(X,1)*size(X,2) % 0 means less data in this pid
            disp(['meta id ' num2str(metaid) ' named ' X.label{3}(metaid,:) ' contains ' num2str(data_count) ' data, out of ' num2str(size(X,1)*size(X,2)) '.'])
            kept_metabolite(metaid) = 0;
            do_loop = 1;
        end
    end
    
    kept_metabolite_id = metabolite_id(find(kept_metabolite==1));
    X = X(:,:,kept_metabolite_id);
end
end


