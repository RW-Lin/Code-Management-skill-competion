function [PAT,matched_PPG_peaks,peaks_ecg_time] = calculate_PAT(peaks_ecg_time,peaks_ppg_time,psqi_session,fs_ppg,timestamp,lag_before,sbp,debug)
    if nargin == 6, debug = 0; sbp = [];end
    
    %% 
    % This script calculates pulse arrival time
    % peaks_ecg_time: in seconds
    % peaks_ppg_time: in seconds
    % psqi_session: signal quanlity index for each pulse peak
    % timestamp: for supine stand test in seconds
    % lag_before: lag before test

    % Author: Runwei Lin
    % Date: 06-12-2024

   %% 
   psqi_session = psqi_session(:);
   matched_PPG_peaks = NaN(size(peaks_ecg_time));  matched_PPG_idx = [];
   unique_ppg_peaks = [];  unique_ecg_peaks = []; unique_matched_PPG_idx = [];

    % here we match pulse peaks to each R peak
   for i = 1:length(peaks_ecg_time)
       time_diff = (peaks_ppg_time(:) - peaks_ecg_time(i));
       
       valid_indices = find(time_diff>0 & time_diff <= 0.5);
      
       if ~isempty(valid_indices)
           [~, min_idx] = min(time_diff(valid_indices));
           matched_PPG_peaks(i) = peaks_ppg_time(valid_indices(min_idx));
           matched_PPG_idx = [matched_PPG_idx [valid_indices(min_idx);i]];
       end
   end


    % remove duplicates: each R peaks corresopnds to one Pulse peak
    
    % PAT
    PAT = 1e3*(matched_PPG_peaks - peaks_ecg_time); % [ms]
    HR = 1./diff(peaks_ecg_time); PR = 1./diff(peaks_ppg_time);

    matched_PPG_peaks = matched_PPG_peaks./60;
    peaks_ecg_time = peaks_ecg_time./60;
    if debug
        figure, subplot(311), plot(matched_PPG_peaks,PAT,'*-')
        ylabel('PAT [ms]'); xline(timestamp,'r')
        subplot(312), plot(matched_PPG_peaks(1:end-1),diff(PAT),'*-')
        ylabel('PATv [ms]'); xline(timestamp,'r');
        subplot(313)
        plot(peaks_ecg_time(1:end-1),HR,'o-','LineWidth',1.5);hold on
        plot(peaks_ppg_time(1:end-1),PR,'*-')
        plot(matched_PPG_peaks,max(HR)*PAT./max(PAT),'*-'); legend('HR','PR','PAT')
        xline(timestamp,'r');
    end
end
