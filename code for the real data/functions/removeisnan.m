function X = removeisnan(X)

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