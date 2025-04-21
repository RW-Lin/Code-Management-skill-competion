function ranksum_pvalues = consistent_analysis_age(feature_old,feature_young,debug)
if nargin == 2, debug = 0; end
% close all
name_fields = fieldnames(feature_old);
ranksum_pvalues = zeros(1,length(name_fields)); 
%% outputs
% figure,
% for nn = 1:length(name_fields)
%    subplot(5,4,nn);plot(feature_young.(name_fields{nn}),'o-');hold on;plot(feature_old.(name_fields{nn}),'*-');
%    title(name_fields{nn});
% end 
% legend('young','old')




%% ranksum
for nn = 1:length(name_fields)
    if (swtest(feature_young.(name_fields{nn})) || swtest(feature_old.(name_fields{nn})))
        [p_nn,h] = ranksum(feature_young.(name_fields{nn}),feature_old.(name_fields{nn}));
        ranksum_pvalues(nn) = p_nn;
    else
        [~,p_nn] = ttest2(feature_young.(name_fields{nn}),feature_old.(name_fields{nn}));
        ranksum_pvalues(nn) = p_nn;
    end
end 

end