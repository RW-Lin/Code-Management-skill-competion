function plot_figures_hrvprv(data_HRV, data_PRV, group)
%% 数据加载
if strcmp(group, 'young')
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

%% 子图设置
subfieldnames = fieldnames(PRV_supine);
phasename = {'supine', 'transition' ,'upright'};

figure;
tiledlayout(ceil(length(subfieldnames) / 5), 5);

% 循环绘制每个指标的箱线图
for subf = 1:length(subfieldnames)
    name = subfieldnames{subf};
    nexttile;
    Compare_boxplot(eval(['HRV_' phasename{1} '.' name]), eval(['HRV_' phasename{2} '.' name]),...
        eval(['HRV_' phasename{3} '.' name]), eval(['PRV_' phasename{1} '.' name]),...
        eval(['PRV_' phasename{2} '.' name]), eval(['PRV_' phasename{3} '.' name]), name);
end

% 添加通用标签


end

function Compare_boxplot(hrv_p1, hrv_p2, hrv_p3, prv_p1, prv_p2, prv_p3, name)
    % 数据整理
    hrv_data = {hrv_p1(:), hrv_p2(:), hrv_p3(:)};
    prv_data = {prv_p1(:), prv_p2(:), prv_p3(:)};
    
    % 确定最大数据长度并填充 NaN
    max_length = max(cellfun(@length, [hrv_data, prv_data]));
    hrv_data_padded = cellfun(@(x) [x; nan(max_length - length(x), 1)], hrv_data, 'UniformOutput', false);
    prv_data_padded = cellfun(@(x) [x; nan(max_length - length(x), 1)], prv_data, 'UniformOutput', false);
    nan_data = nan(size(horzcat(prv_data_padded{:})));
    
    hrv_matrix = horzcat(hrv_data_padded{:});
    prv_matrix = horzcat(prv_data_padded{:});
    % 将数据展开为长格式
    data_all = [horzcat(hrv_data_padded{:}), horzcat(prv_data_padded{:}), nan_data];
    data_all = data_all(:); % 转换为列向量

    % 创建分组标签
    num_phases = 3; % Supine, Transition, Upright
    num_types = 3; % HRV, PRV
    total_points = numel(data_all);

    % 为每个数据点创建对应的分类标签
    phase_labels = repmat({'Supine', 'Transition', 'Upright'}, max_length, num_types);
    phase_labels = phase_labels(:); % 展开为列向量
    
    type_labels = [repmat({'HRV'}, max_length * num_phases, 1); ...
                   repmat({'PRV'}, max_length * num_phases, 1); ...
                   repmat({'1'}, max_length * num_phases, 1)];

    % 转换为 categorical 类型
    phase_labels = categorical(phase_labels, {'Supine', 'Transition', 'Upright'});
    type_labels = categorical(type_labels, {'HRV', 'PRV','1'});

    % 绘制箱线图
    % figure;
    hold on;
    
    % 配对折线部分（连接每阶段的 HRV 和 PRV 数据）
    num_subjects = size(hrv_matrix, 1);
    num_phases = size(hrv_matrix, 2);
    for phase = 1:2:5
        for i = 1:num_subjects
            % 用虚线连接同一阶段的 HRV 和 PRV 数据
            plot([phase - 0.45, phase + 0.45], [hrv_matrix(i, (phase+1)/2), prv_matrix(i, (phase+1)/2)], ...
                 'Color', [0.6, 0.6, 0.6, 0.5], 'LineWidth', 1.2); % 浅灰色线
            % 散点标注
            scatter(phase - 0.45, hrv_matrix(i, (phase+1)/2), 60, 'b', 'filled'); % HRV 点（蓝色）
            scatter(phase + 0.45, prv_matrix(i, (phase+1)/2), 60, 'r', 'filled'); % PRV 点（红色）
        end
    end
    
    % 箱线图部分
    for phase = 1:2:5
        % HRV 的箱线图
        boxchart((phase - 0.45) * ones(num_subjects, 1), hrv_matrix(:, (phase+1)/2), ...
                 'BoxFaceColor', [0, 0, 1], 'LineWidth', 1.5); % 蓝色箱线图
        % PRV 的箱线图
        boxchart((phase + 0.45) * ones(num_subjects, 1), prv_matrix(:, (phase+1)/2), ...
                 'BoxFaceColor', [1, 0, 0], 'LineWidth', 1.5); % 红色箱线图
    end

    % 设置图形属性
    xticks(1:5);
    xticklabels({'Supine', '', 'Transition','', 'Upright'}); % 阶段标签
    
    ylabel(name);
    if strcmp(name,'SD1')   
        legend({'','HRV', 'PRV'}, 'Location', 'Best'); legend boxoff
    end
end
