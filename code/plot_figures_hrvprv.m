function plot_figures_hrvprv(data_HRV,data_PRV,group)
%%
if strcmp(group,'young')
HRV_transition = data_HRV.HRV_transition_young; HRV_supine = data_HRV.HRV_supine_young;
HRV_upright = data_HRV.HRV_upright_young;

PRV_upright = data_PRV.PRV_upright_young; PRV_transition = data_PRV.PRV_transition_young;
PRV_supine = data_PRV.PRV_supine_young; 
else
HRV_transition = data_HRV.HRV_transition_old; HRV_supine = data_HRV.HRV_supine_old;
HRV_upright = data_HRV.HRV_upright_old;

PRV_upright = data_PRV.PRV_upright_old; PRV_transition = data_PRV.PRV_transition_old;
PRV_supine = data_PRV.PRV_supine_old; 
end
%% boxplot 

subfieldnames = fieldnames(PRV_supine); phasename = {'supine','transition','upright'};

figure, t = tiledlayout(3, 5);
for subf = 1:length(subfieldnames)

name = subfieldnames{subf};
Compare_boxplot(eval(['HRV_' phasename{1} '.' name]), eval(['HRV_' phasename{2} '.' name]),...
    eval(['HRV_' phasename{3} '.' name]), eval(['PRV_' phasename{1} '.' name]),...
    eval(['PRV_' phasename{2} '.' name]), eval(['PRV_' phasename{3} '.' name]),name);
if subf>10, xlabel('HRV          PRV','FontSize',10);end
end

end
function Compare_boxplot(hrv_p1,hrv_p2,hrv_p3, prv_p1,prv_p2,prv_p3,name)
    
    hrv_p1 = hrv_p1(:); hrv_p2 = hrv_p2(:);
    hrv_p3 = hrv_p3(:); prv_p1 = prv_p1(:);
    prv_p2 = prv_p2(:); prv_p3 = prv_p3(:);
    
    max_length = max([length(hrv_p1), length(hrv_p2),...
        length(hrv_p3),length(prv_p1), length(prv_p2),...
        length(prv_p3)]);
    hrv_p1_padded = [hrv_p1; nan(max_length - length(hrv_p1), 1)];
    hrv_p2_padded = [hrv_p2; nan(max_length - length(hrv_p2), 1)];
    hrv_p3_padded = [hrv_p3; nan(max_length - length(hrv_p3), 1)];
    prv_p1_padded = [prv_p1; nan(max_length - length(prv_p1), 1)];
    prv_p2_padded = [prv_p2; nan(max_length - length(prv_p2), 1)];
    prv_p3_padded = [prv_p3; nan(max_length - length(prv_p3), 1)];
    Nanpadd = nan(size(prv_p3_padded));
    % 合并为矩阵
    data_matrix = [hrv_p1_padded, hrv_p2_padded,hrv_p3_padded,Nanpadd...
        prv_p1_padded, prv_p2_padded,prv_p3_padded];
    hrv_matrix = [hrv_p1_padded, hrv_p2_padded,hrv_p3_padded];
    prv_matrix = [prv_p1_padded, prv_p2_padded,prv_p3_padded];
    nexttile;
    boxplot(data_matrix, 'Labels', {'Supine', 'Transition','Standing','','Supine', 'Transition','Standing'},'Notch','on');
    ax = gca;
    ax.FontSize=7;
    ax.XTickLabelRotation = 60; % 调整标签角度
    ylabel(name,'FontSize',10);

    %% statistical analysis + correction -> label
    % 
    % [p,~,stats] = kruskalwallis(hrv_matrix, {'Supine-HRV', 'Transition-HRV','Upright-HRV'});
    % [c,m,h] = multcompare(stats,"CriticalValueType","dunn-sidak");
    % [p,~,stats] = kruskalwallis(hrv_matrix, {'Supine-PRV', 'Transition-PRV','Upright-PRV'});
    % [c,m,h] = multcompare(stats,"CriticalValueType","dunn-sidak");
    
end
