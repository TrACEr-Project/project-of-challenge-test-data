%% This is for choosing data type: T0 , T0corrected or full dynamic data?

function X=generate_data_types(X, data_type)


switch data_type
    case 'T0'
        X=X(:,1,:);
    case 'T0 corrected'
        X=X(:,2:end,:)-X(:,1,:);
    case 'full dynamic'
        X=X;
end
end







