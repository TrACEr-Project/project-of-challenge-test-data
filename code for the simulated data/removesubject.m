
%%%%%% This code is to remove the outlier in the first (subject) mode 

function Y=removeoutlier(Y,outlier_pid)

[subject_notoutlier,index] = setdiff(str2num(Y.label{1}), outlier_pid);
Y=Y(index,:,:);
end


