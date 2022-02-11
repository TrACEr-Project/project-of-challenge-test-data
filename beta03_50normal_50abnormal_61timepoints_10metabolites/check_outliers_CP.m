%%% This function is about checking the outliers for the CP model
function check_outliers_CP(X, Fac)
W=ones(size(X));
W(find(isnan(X.data)))=0;
%% residual
residual_temp = (tensor(X) - full(Fac)).*W;
residual_now = sum(residual_temp.data.^2,[2,3]);
% residual_now = residual_now/sum(residual_now);
% residual_std = std(residual_now);
% leverage
A = Fac.U{1};
leverage_now = diag(A*inv(A'*A)*A');

figure
plot(leverage_now,residual_now,'+')
xlabel('Leverage score')
ylabel('Residual')
title('subjects')

for i=1:size(X,1)
    Ind_num{i} = num2str(i);
end
text(leverage_now,residual_now,Ind_num)



%% residual in metabolites
residual_temp_meta = (tensor(X) - full(Fac)).*W;
residual_now_meta = sum(residual_temp.data.^2,[1,3]);
% residual_now_meta = residual_now_meta/sum(residual_now_meta);

% leverage
B = Fac.U{2};
leverage_now_meta = diag(B*inv(B'*B)*B');
labelss = {'Ins_B','GLC_B','Pyr_B','Lac_B','Ala_B','Glyc_B','Bhb_B','FFA_B',...
           'TG_B','Chol_B'};
figure
plot(leverage_now_meta,residual_now_meta,'bo','MarkerFaceColor','b')
xlabel('Leverage score')
ylabel('Residual')
title('metabolites')

% for i=1:size(X,2)
%     Ind_num_meta{i} = num2str(i);
% end

text(leverage_now_meta,residual_now_meta,labelss,'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize', 13)
set(gca,'Fontsize',13)


