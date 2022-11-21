% This script shows how to numerically check the uniqueness of the CP factorization

%%%%%% Input
% Fac_X: the factors obtained using different initializations
% erF: the error of the function value for each initialization
% goodness_X1: the convergence information
%%%%%% Output
% unique_index=0 -> NOT unique
% unique_index=1 -> Unique
% unique_index=2 -> Inconclusive, need more random starts

function unique_index=unique_test_CP(Fac_X, erF, goodness_X1)

switch nargin
    case 3
        good_flag = find(goodness_X1(:) == 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL.' | goodness_X1(:) == 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH.');
    case 2
        good_flag =1: length(erF)
    otherwise
        display('Error--input information')
        return
end

if length(good_flag)>=1
    F_round = round(erF(good_flag),7);
    best_F_index = good_flag(F_round == min(F_round));
    if length(best_F_index) < 2
        F_round = round(erF(good_flag),5);
        best_F_index = good_flag(F_round == min(F_round));
    end
else
    F_round = round(goodness_X(:,2),10);
    best_F_index = find(F_round == min(F_round));
end

eps = .05;
if length(best_F_index)==1
    unique_index = 2;
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
    if worst_check < (1-eps) 
        unique_index = 0;
    else
        unique_index = 1;
    end
end
end