
function NMR=select_metabolites_simul(NMR)

a=NMR.label{3,1};
b=['Insulin        ';'Glucose        ';'Pyruvate       ';'Lactate        ';'Ala            ';'Glycerol       ';'bOHbutyrate    '];
% do not include Glycerol since most of them are missing. Missing check
% will remove this metabolite
selected_index=find( ismember(a, b, 'rows'));
NMR=NMR(:,:,selected_index);
end

