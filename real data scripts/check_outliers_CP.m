%%% This function is about checking the outliers for the CP model
function check_outliers_CP(X, Fac)

%% residual
X_predic=full(Fac);
residual_temp = X.data - X_predic.data;
residual_now =squeeze( nansum(residual_temp.^2,[2,3]));
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
residual_now_meta =squeeze(nansum(residual_temp.^2,[1,2]));


% residual_now_meta = residual_now_meta/sum(residual_now_meta);

% leverage
B = Fac.U{3};
leverage_now_meta = diag(B*inv(B'*B)*B');
labelss = {'Ins','GLC','Pyr','Lac','Ala','Bhb'};
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


