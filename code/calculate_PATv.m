function [PAT,peaks_ecg_time,matched_PPG_peaks] = calculate_PATv(peaks_ecg_time,peaks_ppg_time,psqi_session,fs_ppg,timestamp,lag_before)
    %% 
    % This script calculates pulse arrival time
    % peaks_ecg_time: in seconds
    % peaks_ppg_time: in seconds
    % psqi_session: signal quanlity index for each pulse peak
    % timestamp: for supine stand test
    % lag_before: lag before test

    % Author: Runwei Lin
    % Date: 06-12-2024
   if nargin == 6, debug = 0;end
   %% 
   psqi_session = psqi_session(:);
   matched_PPG_peaks = NaN(size(peaks_ecg_time));  matched_PPG_idx = [];
   unique_ppg_peaks = [];  unique_ecg_peaks = []; unique_matched_PPG_idx = [];

    % here we match pulse peaks to each R peak
   for i = 1:length(peaks_ecg_time)
       time_diff = (peaks_ppg_time(:) - peaks_ecg_time(i));
       
       valid_indices = find(time_diff>0.05 & time_diff <= 0.5);
      
       if ~isempty(valid_indices)
           [~, min_idx] = min(time_diff(valid_indices));
           matched_PPG_peaks(i) = peaks_ppg_time(valid_indices(min_idx));
           matched_PPG_idx = [matched_PPG_idx ;valid_indices(min_idx), i];
       end
   end
    
    %% % unmatched R-peak / PPG peak
    % checked_ppg_peaks = matched_PPG_peaks(~isnan(matched_PPG_peaks));
    % checked_ecg_peaks = peaks_ecg_time(~isnan(matched_PPG_peaks));
    % matched_PPG_idx = matched_PPG_idx(~isnan(matched_PPG_peaks)); %

    % unmatched_ecg_peaks = setdiff(1:length(peaks_ecg_time),matched_PPG_idx);
    % 
    % % duplicate PPG peaks
    % [~, first_idx] = unique(matched_PPG_idx(1,:), 'first'); % 
    % duplicates = setdiff(1:length(matched_PPG_idx), first_idx); %     
    % % 
    % 
    % unique_ppg_indices = unique(matched_PPG_idx(1, :), 'stable');  % 保留首次匹配的PPG峰
    % duplicates = setdiff(matched_PPG_idx(:, 1), unique_ppg_indices);  % 找出重复的PPG峰
    % 
    % unmatched_ecg_idx = find(isnan(matched_PPG_peaks));  % 找出未匹配的R峰
    % % remove duplicates
    % while ~isempty(duplicates)
    %     % find duplicate peaks
    %     duplicate_ppg_peak = matched_PPG_idx(1,duplicates(1));
    % 
    %     % 
    %     duplicate_indices = find(matched_PPG_idx(1,:) == duplicate_ppg_peak);
    % 
    %     % selection
    %     min_time_diff = inf;
    %     best_index = -1;
    %     for idx = duplicate_indices'
    %         time_diff = abs(matched_PPG_idx(1,idx) - checked_ecg_peaks(idx));
    %         if time_diff < min_time_diff
    %             min_time_diff = time_diff;
    %             best_index = idx;
    %         end
    %     end
    % 
    %     % 将时间差最小的 PPG 峰及其对应的 ECG 峰加入去重结果
    %     unique_ppg_peaks = [unique_ppg_peaks; unique(checked_ppg_peaks(best_index))];
    %     unique_ecg_peaks = [unique_ecg_peaks; unique(checked_ecg_peaks(best_index))];
    %     unique_matched_PPG_idx = [unique_matched_PPG_idx; matched_PPG_idx(best_index)];
    %     % 从 duplicates 和 checked_ppg_peaks 中移除该重复 PPG 峰的所有索引
    %     duplicates = setdiff(duplicates, duplicate_indices);
    %     checked_ppg_peaks(duplicate_indices) = NaN; % 标记为 NaN
    % end
    % 
    % % 将剩余唯一的非重复 PPG 峰和 ECG 峰加入结果
    % remaining_indices = find(~isnan(checked_ppg_peaks));
    % unique_ppg_peaks = sort([unique_ppg_peaks checked_ppg_peaks(remaining_indices)]);
    % unique_ecg_peaks = sort([unique_ecg_peaks checked_ecg_peaks(remaining_indices)]);
    
    %% Remove duplicates
    unique_ppg_indices = unique(matched_PPG_idx(:, 1), 'stable');  % 保留首次匹配的PPG峰
    duplicates = setdiff(matched_PPG_idx(:,1), unique_ppg_indices);  % 找出重复的PPG峰
    unmatched_ecg_idx = setdiff(1:length(peaks_ecg_time), matched_PPG_idx(:,2));
    
    for duplicate = duplicates'
        % duplicate_indices = find(matched_PPG_idx(:, 1) == duplicate);

        % 
        best_r_idx = -1;
        min_time_diff = inf;

        for ecg_idx = unmatched_ecg_idx'
            time_diff = abs(peaks_ppg_time(duplicate) - peaks_ecg_time(ecg_idx));
            if time_diff < min_time_diff && time_diff <= 0.5
                min_time_diff = time_diff;
                best_r_idx = ecg_idx;
            end
        end

        % Update
        if best_r_idx ~= -1
            matched_PPG_peaks(best_r_idx) = peaks_ppg_time(duplicate);
            matched_PPG_idx = [matched_PPG_idx; duplicate, best_r_idx];
            unmatched_ecg_idx(unmatched_ecg_idx == best_r_idx) = [];  % 移除已匹配的R峰
        end
    end

    %% remove the scenario one ppg peaks accounting for several rpeak
    [unique_ppg_indices, ~, ic] = unique(matched_PPG_idx(:, 1), 'stable');  % 找出PPG峰索引
    counts = accumarray(ic, 1);  % 统计每个PPG峰的出现次数
    duplicates = unique_ppg_indices(counts > 1);  % 重复的PPG峰索引
    
    % 初始化未匹配的R峰列表
    unmatched_rpeaks = 1:length(peaks_ecg_time);

    % 遍历重复的PPG峰
    for i = 1:length(duplicates)
        duplicate_ppg = duplicates(i);  % 当前重复的PPG峰
        rows = find(matched_PPG_idx(:, 1) == duplicate_ppg);  % 所有匹配到该PPG峰的R峰行
        
        % 计算时间差，选择时间差最小的R峰作为保留项
        time_diffs = abs(peaks_ppg_time(duplicate_ppg) - peaks_ecg_time(matched_PPG_idx(rows, 2)));
        time_diffs(time_diffs>0.5 | time_diffs<0.05) = Inf;
        [~, min_idx] = min(time_diffs);
        best_row = rows(min_idx);
        
        % 保留最佳匹配
        best_rpeak_idx = matched_PPG_idx(best_row, 2);
        matched_PPG_peaks(best_rpeak_idx) = peaks_ppg_time(duplicate_ppg);
        unmatched_rpeaks(unmatched_rpeaks == best_rpeak_idx) = [];  % 移除已匹配的R峰
    
        % 处理其他R峰：重新分配剩余的PPG峰
        rows(rows == best_row) = [];  % 剩余需要重新分配的R峰
        for j = 1:length(rows)
            rpeak_idx = matched_PPG_idx(rows(j), 2);
            
            % 在未匹配的PPG峰中查找时间差最小的
            remaining_ppg_indices = setdiff(1:length(peaks_ppg_time), matched_PPG_idx(:, 1), 'stable');
            time_diffs = abs(peaks_ppg_time(remaining_ppg_indices) - peaks_ecg_time(rpeak_idx));
            valid_indices = find(time_diffs > 0.05 & time_diffs <= 0.5);
            
            if ~isempty(valid_indices)
                % 选择时间差最小的PPG峰
                [~, min_idx] = min(time_diffs(valid_indices));
                best_ppg_idx = remaining_ppg_indices(valid_indices(min_idx));
                matched_PPG_peaks(rpeak_idx) = peaks_ppg_time(best_ppg_idx);
                matched_PPG_idx = [matched_PPG_idx; best_ppg_idx, rpeak_idx];
            else
                % 无法重新匹配，设为NaN
                matched_PPG_peaks(rpeak_idx) = NaN;
            end
        end
    end
    %% 
    % PAT
    PAT = 1e3*(matched_PPG_peaks - peaks_ecg_time);  
    if debug,
        figure, subplot(411), plot(matched_PPG_peaks,PAT,'*-')
        ylabel('PAT [ms]'); xline(timestamp,'r')
        subplot(412), plot(matched_PPG_peaks(1:end-1),diff(PAT),'*-')
        ylabel('PATv [ms]'); xline(timestamp,'r');
        subplot(413)
        plot(peaks_ecg_time(1:end-1),60./diff(peaks_ecg_time),'o-','LineWidth',1.5);hold on
        plot(peaks_ppg_time(1:end-1),60./diff(peaks_ppg_time),'*-');
        plot(matched_PPG_peaks,PAT,'*-'); legend('HR','PR','PAT');  xline(timestamp,'r');
        subplot(414)
        plot(peaks_ecg_time(1:end-1),diff(peaks_ecg_time),'o-','LineWidth',1.5);hold on
        plot(peaks_ppg_time(1:end-1),diff(peaks_ppg_time),'*-');hold off
        ylabel('R-R [ms]');xline(timestamp,'r');
    end
    matched_PPG_peaks = matched_PPG_peaks./60; peaks_ecg_time = peaks_ecg_time./60;
    %% spectral analysis PAT
    % Fs_resample = 4;
    % new_PAT_time = min(matched_PPG_peaks)*60:1/Fs_resample:max(matched_PPG_peaks)*60;
    % new_PAT = interp1(60*matched_PPG_peaks(~isnan(PAT)), PAT(~isnan(PAT)), new_PAT_time,'pchip');
    % if ~debug
    %     figure,subplot(211)
    %     plot(matched_PPG_peaks*60,PAT,'*-');hold on
    %     plot(new_PAT_time,new_PAT,'s-')
    %     subplot(212)
    %     wvd(new_PAT,Fs_resample,'smoothedPseudo')
    % end
end
