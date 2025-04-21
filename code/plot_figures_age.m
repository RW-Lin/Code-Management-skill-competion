function plot_figures_age(data_HRV,type,activity_type)
% plots for EMBC 2025
%%
if strcmp(type,'HRV')
    HRV_transition_old = data_HRV.HRV_transition_old; HRV_transition_young = data_HRV.HRV_transition_young;
    HRV_supine_old = data_HRV.HRV_supine_old; HRV_supine_young = data_HRV.HRV_supine_young;
    HRV_upright_old = data_HRV.HRV_upright_old;  HRV_upright_young = data_HRV.HRV_upright_young;
else
    HRV_transition_old = data_HRV.PRV_transition_old; HRV_transition_young = data_HRV.PRV_transition_young;
    HRV_supine_old = data_HRV.PRV_supine_old; HRV_supine_young = data_HRV.PRV_supine_young;
    HRV_upright_old = data_HRV.PRV_upright_old;  HRV_upright_young = data_HRV.PRV_upright_young;
end
%%
phasename = {'supine','transition','upright'};

if strcmp(type,'HRV')
subfieldnames = fieldnames(HRV_transition_old); 
figure, t = tiledlayout(5, 5);
for subf = 1:length(subfieldnames)
name = subfieldnames{subf};
Compare_boxplot(eval(['HRV_' phasename{1} '_old.' name]), eval(['HRV_' phasename{2} '_old.' name]),...
    eval(['HRV_' phasename{3} '_old.' name]), eval(['HRV_' phasename{1} '_young.' name]),...
    eval(['HRV_' phasename{2} '_young.' name]), eval(['HRV_' phasename{3} '_young.' name]),name);
if subf>15, xlabel('Old          Young','FontSize',11);end
end
% figure, t =  tiledlayout(2, 2);
% for subf = 1:length(subfieldnames)-4
% name = subfieldnames{subf};
% Compare_boxplot(eval(['HRV_' phasename{1} '_old.' name]), eval(['HRV_' phasename{2} '_old.' name]),...
%     eval(['HRV_' phasename{3} '_old.' name]), eval(['HRV_' phasename{1} '_young.' name]),...
%     eval(['HRV_' phasename{2} '_young.' name]), eval(['HRV_' phasename{3} '_young.' name]),name);
% if subf>10, xlabel('Old          Young','FontSize',11);end
% end
else
subfieldnames = fieldnames(HRV_transition_old); 
figure, t = tiledlayout(4, 5);
for subf = 1:length(subfieldnames)
name = subfieldnames{subf};
Compare_boxplot(eval(['HRV_' phasename{1} '_old.' name]), eval(['HRV_' phasename{2} '_old.' name]),...
    eval(['HRV_' phasename{3} '_old.' name]), eval(['HRV_' phasename{1} '_young.' name]),...
    eval(['HRV_' phasename{2} '_young.' name]), eval(['HRV_' phasename{3} '_young.' name]),name);
if subf>10, xlabel('Old          Young','FontSize',11);end    
end
end
sgtitle([type '-' activity_type])

%% boxplot time features
% freq domain
% figure,
% subplot(2,3,1)
% Compare_boxplot(HRV_supine_young.LF, HRV_transition_young.LF, HRV_upright_young.LF,...
%     HRV_supine_old.LF, HRV_transition_old.LF, HRV_upright_old.LF,'LF',type);
% subplot(2,3,2)
% Compare_boxplot(HRV_supine_young.HF, HRV_transition_young.HF, HRV_upright_young.HF, ...
%     HRV_supine_old.HF, HRV_transition_old.HF, HRV_upright_old.HF,'HF',type)
% subplot(2,3,3)
% Compare_boxplot(HRV_supine_young.LFHFratio, HRV_transition_young.LFHFratio, HRV_upright_young.LFHFratio, ...
%     HRV_supine_old.LFHFratio, HRV_transition_old.LFHFratio, HRV_upright_old.LFHFratio,'LFHF',type)
% subplot(2,3,4)
% Compare_boxplot(HRV_supine_young.pLF, HRV_transition_young.pLF, HRV_upright_young.pLF, ...
%     HRV_supine_old.pLF, HRV_transition_old.pLF, HRV_upright_old.pLF,'pLF',type)
% subplot(2,3,5)
% Compare_boxplot(HRV_supine_young.pHF, HRV_transition_young.pHF, HRV_upright_young.pHF, ...
%     HRV_supine_old.pHF, HRV_transition_old.pHF, HRV_upright_old.pHF,'pHF',type)
% subplot(2,3,6)
% Compare_boxplot(HRV_supine_young.SDratio, HRV_transition_young.SDratio, HRV_upright_young.SDratio,...
%     HRV_supine_old.SDratio, HRV_transition_old.SDratio, HRV_upright_old.SDratio,'SDratio',type)
% sgtitle(activity_type)

% Time domain
% figure,
% subplot(2,3,1)
% Compare_boxplot(HRV_supine_young.SDNN, HRV_transition_young.SDNN, HRV_upright_young.SDNN,...
%     HRV_supine_old.SDNN, HRV_transition_old.SDNN, HRV_upright_old.SDNN,'SDNN',type)
% subplot(2,3,2)
% Compare_boxplot(HRV_supine_young.SDSD, HRV_transition_young.SDSD, HRV_upright_young.SDSD, ...
%     HRV_supine_old.SDSD, HRV_transition_old.SDSD, HRV_upright_old.SDSD,'SDSD',type)
% subplot(2,3,3)
% Compare_boxplot(HRV_supine_young.RMSSD, HRV_transition_young.RMSSD, HRV_upright_young.RMSSD, ...
%     HRV_supine_old.RMSSD, HRV_transition_old.RMSSD, HRV_upright_old.RMSSD,'RMSSD',type)
% subplot(2,3,4)
% Compare_boxplot(HRV_supine_young.pNN50, HRV_transition_young.pNN50, HRV_upright_young.pNN50, ...
%     HRV_supine_old.pNN50, HRV_transition_old.pNN50, HRV_upright_old.pNN50,'pNN50',type)
% subplot(2,3,5)
% Compare_boxplot(HRV_supine_young.SD1, HRV_transition_young.SD1, HRV_upright_young.SD1, ...
%     HRV_supine_old.SD1, HRV_transition_old.SD1, HRV_upright_old.SD1,'SD1',type)
% subplot(2,3,6)
% Compare_boxplot(HRV_supine_young.SD2, HRV_transition_young.SD2, HRV_upright_young.SD2,...
%     HRV_supine_old.SD2, HRV_transition_old.SD2, HRV_upright_old.SD2,'SD2',type)
% sgtitle(activity_type)

% freq domain

% PAT
% if isfield(HRV_transition_old,"PATmn")
%     figure,
%     subplot(2,3,1)
%     Compare_boxplot(HRV_supine_young.PATmn, HRV_transition_young.PATmn, HRV_upright_young.PATmn,...
%         HRV_supine_old.PATmn, HRV_transition_old.PATmn, HRV_upright_old.PATmn,'PATmn',type);
%     subplot(2,3,2)
%     Compare_boxplot(HRV_supine_young.PATsd, HRV_transition_young.PATsd, HRV_upright_young.PATsd, ...
%         HRV_supine_old.PATsd, HRV_transition_old.PATsd, HRV_upright_old.PATsd,'PATsd',type)
%     subplot(2,3,3)
%     Compare_boxplot(HRV_supine_young.PATvmn, HRV_transition_young.PATvmn, HRV_upright_young.PATvmn, ...
%         HRV_supine_old.PATvmn, HRV_transition_old.PATvmn, HRV_upright_old.PATvmn,'PATvmn',type)
%     subplot(2,3,4)
%     Compare_boxplot(HRV_supine_young.PATvsd, HRV_transition_young.PATvsd, HRV_upright_young.PATvsd, ...
%         HRV_supine_old.PATvsd, HRV_transition_old.PATvsd, HRV_upright_old.PATvsd,'PATvsd',type)
%     subplot(2,3,5)
%     Compare_boxplot(HRV_supine_young.AVNN, HRV_transition_young.AVNN, HRV_upright_young.AVNN, ...
%         HRV_supine_old.AVNN, HRV_transition_old.AVNN, HRV_upright_old.AVNN,'AHR',type)
% 
% end

end
function Compare_boxplot(data_p1_young,data_p2_young,data_p3_young,...
    data_p1_old,data_p2_old,data_p3_old,name,type)
    
    data_p1_young = data_p1_young(:); data_p2_young = data_p2_young(:);
    data_p3_young = data_p3_young(:); data_p1_old = data_p1_old(:);
    data_p2_old = data_p2_old(:); data_p3_old = data_p3_old(:);
    
    max_length = max([length(data_p1_young), length(data_p2_young),...
        length(data_p3_young),length(data_p1_old), length(data_p2_old),...
        length(data_p3_old)]);
    data_p1_young_padded = [data_p1_young; nan(max_length - length(data_p1_young), 1)];
    data_p2_young_padded = [data_p2_young; nan(max_length - length(data_p2_young), 1)];
    data_p3_young_padded = [data_p3_young; nan(max_length - length(data_p3_young), 1)];
    data_p1_old_padded = [data_p1_old; nan(max_length - length(data_p1_old), 1)];
    data_p2_old_padded = [data_p2_old; nan(max_length - length(data_p2_old), 1)];
    data_p3_old_padded = [data_p3_old; nan(max_length - length(data_p3_old), 1)];
    
    % 合并为矩阵
    data_matrix = [data_p1_young_padded, data_p2_young_padded,data_p3_young_padded,...
        data_p1_old_padded, data_p2_old_padded,data_p3_old_padded];
    nexttile;
    boxplot(data_matrix, 'Labels', {'Supine', 'Transition','Upright','Supine', 'Transition','Upright'},'FactorSeparator', 1);
    ax = gca;
    ax.FontSize=7;
    ax.XTickLabelRotation = 60; % 调整标签角度
    ylabel(name);

    %% statistical analysis + correction -> label
    % 
    % [p,~,stats] = kruskalwallis(data_matrix, {'Supine-Young', 'Transition-Young','Upright-Young',...
    %     'Supine-Old', 'Transition-Old','Upright-Old'});
    % [c,m,h] = multcompare(stats,"CriticalValueType","dunn-sidak");

    
end
