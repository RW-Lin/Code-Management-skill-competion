function output = consistent_analysis(HRV_feature,PRV_feature,debug)
if nargin == 2, debug = 0; end
% close all
name_fields = fieldnames(PRV_feature);
pvalues = zeros(1,length(name_fields)); spcc = zeros(1,length(name_fields)); 
%% outputs
% figure,
% for nn = 1:length(name_fields)
%    subplot(5,4,nn);plot(HRV_feature.(name_fields{nn}),'o-');hold on;plot(PRV_feature.(name_fields{nn}),'*-');
%    title(name_fields{nn});
% end 
% legend('hrv','prv')
%% ranksum
for nn = 1:length(name_fields)
    
    % plot_HRV_PRV_correlation(HRV_feature.(name_fields{nn}), PRV_feature.(name_fields{nn}), name_fields{nn})
    if (swtest(HRV_feature.(name_fields{nn})) && swtest(PRV_feature.(name_fields{nn})))
        [~,p_nn] = ttest(HRV_feature.(name_fields{nn}),PRV_feature.(name_fields{nn}));
        pvalues(nn) = p_nn;
    else
        [p_nn,h] = signrank(HRV_feature.(name_fields{nn}),PRV_feature.(name_fields{nn}));        
        pvalues(nn) = p_nn;
    end
end 
%% Bland-Altman analysis
% LF 参数的一致性分析
% for nn = 1:length(name_fields)
%     if strcmp(name_fields{nn},'LF') | strcmp(name_fields{nn},'HF') | strcmp(name_fields{nn},'LFHFratio')
%         [rpc, fig] = BlandAltman(HRV_feature.(name_fields{nn}).', PRV_feature.(name_fields{nn}).', ...
%                [name_fields{nn}], 'baInfo', {'RPC%', 'ks', 'LOA', 'CV','SSE','RMSE'});
%     end
% end

%% Correlation
% for nn = 1:length(name_fields)
%     spcc(nn) = corr(HRV_feature.(name_fields{nn}).',PRV_feature.(name_fields{nn}).','type','Spearman');
% end
%%
% 初始化变量
nFeatures = length(name_fields);
ranksum_pvalues = zeros(nFeatures, 1); % p-value 记录
feature_stats = struct();              % 用于记录统计量

% 循环每个特征进行 ranksum 检验，并记录统计量
for nn = 1:nFeatures
    % 提取当前特征数据
    HRV_data = HRV_feature.(name_fields{nn});
    PRV_data = PRV_feature.(name_fields{nn});
    
    % 进行 ranksum 检验
    [p_nn, ~] = ranksum(HRV_data, PRV_data);
    ranksum_pvalues(nn) = p_nn;

    % 计算 HRV 特征的统计量
    Q1_HRV = quantile(HRV_data, 0.25); % Q1
    Q3_HRV = quantile(HRV_data, 0.75); % Q3
    IQR_HRV = Q3_HRV - Q1_HRV;         % IQR
    lower_bound_HRV = Q1_HRV - 1.5 * IQR_HRV; % 正确计算 lower_bound
    upper_bound_HRV = Q3_HRV + 1.5 * IQR_HRV;

    % 寻找实际满足边界条件的数据
    HRV_in_bounds = HRV_data(HRV_data >= lower_bound_HRV & HRV_data <= upper_bound_HRV);
    HRV_min_adj = min(HRV_in_bounds); % Lower adjacent value
    HRV_max_adj = max(HRV_in_bounds); % Upper adjacent value

    % 存储 HRV 统计信息
    feature_stats.(name_fields{nn}).HRV.median = median(HRV_data);
    feature_stats.(name_fields{nn}).HRV.Q1 = Q1_HRV;
    feature_stats.(name_fields{nn}).HRV.Q3 = Q3_HRV;
    feature_stats.(name_fields{nn}).HRV.lower_adj = HRV_min_adj; % 实际下界值
    feature_stats.(name_fields{nn}).HRV.upper_adj = HRV_max_adj; % 实际上界值
    feature_stats.(name_fields{nn}).HRV.bounds = [lower_bound_HRV, upper_bound_HRV];

    % 计算 PRV 特征的统计量
    Q1_PRV = quantile(PRV_data, 0.25);
    Q3_PRV = quantile(PRV_data, 0.75);
    IQR_PRV = Q3_PRV - Q1_PRV;
    lower_bound_PRV = Q1_PRV - 1.5 * IQR_PRV;
    upper_bound_PRV = Q3_PRV + 1.5 * IQR_PRV;

    % 寻找实际满足边界条件的数据
    PRV_in_bounds = PRV_data(PRV_data >= lower_bound_PRV & PRV_data <= upper_bound_PRV);
    PRV_min_adj = min(PRV_in_bounds); % Lower adjacent value
    PRV_max_adj = max(PRV_in_bounds); % Upper adjacent value

    % 存储 PRV 统计信息
    feature_stats.(name_fields{nn}).PRV.median = median(PRV_data);
    feature_stats.(name_fields{nn}).PRV.Q1 = Q1_PRV;
    feature_stats.(name_fields{nn}).PRV.Q3 = Q3_PRV;
    feature_stats.(name_fields{nn}).PRV.lower_adj = PRV_min_adj; % 实际下界值
    feature_stats.(name_fields{nn}).PRV.upper_adj = PRV_max_adj; % 实际上界值
    feature_stats.(name_fields{nn}).PRV.bounds = [lower_bound_PRV, upper_bound_PRV];
end

% for i = 1:nFeatures
%     disp('HRV:') 
%     disp(feature_stats.(name_fields{nn}).HRV)
%     disp('PRV:')
%     disp(feature_stats.(name_fields{nn}).PRV)
% end

%%
% output.spcc = spcc; 
output.pvalues = pvalues;
output.ranksum = ranksum_pvalues; output.feature = feature_stats;

end

function plot_HRV_PRV_correlation(HRV_data, PRV_data, feature_name)
    % 创建图形窗口
    figure;

    % 第一子图：散点图
    subplot(1, 2, 1); % 设置为第一个子图
    hold on;
    scatter(HRV_data, PRV_data, 50, 'b', 'filled'); % 绘制散点图
    x_range = [min([HRV_data; PRV_data]), max([HRV_data; PRV_data])];
    plot(x_range, x_range, 'r--', 'LineWidth', 1.5); % 添加45°参考线
    xlabel('HRV Feature');
    ylabel('PRV Feature');
    title(['Scatter Plot: ', feature_name]);
    legend({'Data Points', 'Perfect Agreement Line'}, 'Location', 'Best');
    grid on;
    hold off;

    % 第二子图：折线箱型图
    subplot(1, 2, 2); % 设置为第二个子图
    hold on;
    % 绘制箱型图
    data_to_plot = [HRV_data', PRV_data'];
    boxplot(data_to_plot, 'Labels', {'HRV', 'PRV'}, 'Widths', 0.6);
    % 添加配对数据的折线
    for i = 1:length(HRV_data)
        plot([1, 2], [HRV_data(i), PRV_data(i)], 'k--'); % 配对数据用虚线连接
    end
    % 添加散点
    scatter(ones(size(HRV_data)), HRV_data, 50, 'b', 'filled'); % HRV 数据散点
    scatter(2 * ones(size(PRV_data)), PRV_data, 50, 'r', 'filled'); % PRV 数据散点

    % 图形样式
    xlabel('Feature Type');
    ylabel('Feature Value');
    title(['Paired Boxplot: ', feature_name]);
    legend({'Pairwise Connection', 'HRV Data', 'PRV Data'}, 'Location', 'Best');
    grid on;
    hold off;
end
