%%  Main program compares old / young under different positions
clear all,close all

Fs_resample = 4; %s^-1
overlap = 1; win_len = 2; %minute
lag_before = 5; lag_after = 3; % duration before and after timestamp
fs_ppg = 1e3; debug = 0; fs_ecg = 3e2;

% please change savepath
savepath = 'C:\Users\LinR\OneDrive - University of Twente\work\project\Spectral analysis of HRV, PTT and BP\Writing\EMBC 2025\arxiv_codes\data\'; 
dirInfo = dir(savepath);
%% initialization
HRV_supine_old = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
HRV_upright_old = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
HRV_transition_old = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
PRV_upright_old = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
PRV_supine_old = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
PRV_transition_old =struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
HRV_supine_young = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
HRV_upright_young = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
HRV_transition_young = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
PRV_upright_young = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
PRV_supine_young = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]);
PRV_transition_young = struct('AHR',[],'RMSSD',[],'SDNN',[],'SDSD',[],'pNN50',[],'LF',[],'HF',[],'pLF',[],...
    'pHF',[],'LFHFratio',[],'SD1',[],'SD2',[],'SDratio',[],'ApEn',[],'SampEn',[]); 

%%
for sg_num = 3:length(dirInfo) %[3:13 16:27 29:31 33:length(dirInfo)] % [15 16 27 36 60
    
    filename = dirInfo(sg_num).name;
    id = sscanf(filename, '%*[^0-9]%d'); id = id(1);
    
    load([savepath filename])
    %
    data_ppg = data_segments.ppg; data_ecg = data_segments.ecg;
    
    ppg_peaktime = data_segments.ppgpeak_time; ecg_peaktime = data_segments.ecgpeak_time;
    session_PPG_time = data_segments.ppgtime; session_ECG_time = data_segments.ecgtime;
    ppg_len = length(data_segments.ppg) / 1e3; timestamp = data_segments.timestamp;
    psqi_session = data_segments.psqi; 
    ecg_peakidx = data_segments.ecg_idx; ppg_peakidx = data_segments.ppg_idx;
    
    [~, calibrated_peaks_ecg, calibrated_peaks_ppg] =...
        calculate_PATv(ecg_peaktime*60,ppg_peaktime*60,psqi_session,fs_ppg,timestamp*60,lag_before);
    
    RR = HRV.RRfilter(diff(calibrated_peaks_ecg),60); 
    PR = HRV.RRfilter(diff(calibrated_peaks_ppg),60);
    time_RR = calibrated_peaks_ecg(1:end-1); time_PR = calibrated_peaks_ppg(1:end-1);
    RR_win = RR(time_RR < timestamp(1)-0.5 & time_RR > session_ECG_time(1)+0.5);
    PR_win = PR(time_PR < timestamp(1)-0.5 & time_PR > session_PPG_time(1)+0.5);

    % ******** visual check if the peak detection is correct *********
    if sg_num == 4
        figure,subplot(311)
        plot(session_PPG_time,data_ppg)
        hold on 
        plot(ppg_peaktime,data_ppg(ppg_peakidx),'o')
        title('PPG [n.u.]');xline(timestamp,'Color','g','LineWidth',2)
        xlim([54.35 54.7])
        subplot(312)
        plot(session_ECG_time,data_ecg)
        hold on 
        plot(ecg_peaktime,data_ecg(ecg_peakidx),'o');
        xline(timestamp,'Color','g','LineWidth',2);ylabel('ECG [mV]');
        xlim([54.35 54.7])      
        subplot(313)
        plot(calibrated_peaks_ppg(1:end-1),diff(calibrated_peaks_ppg),'o-');hold on;
        plot(calibrated_peaks_ecg(1:end-1),diff(calibrated_peaks_ecg),'*-');
        legend('PR','HR');ylabel('Heart rate [bpm]')
        xline(timestamp,'Color','g','LineWidth',2)
        xlim([54.35 54.7])      
    end

    %% supine   
    time_RR_win = time_RR(time_RR < timestamp(1)-.5 & time_RR > session_ECG_time(1)+.5);
    time_PR_win = time_PR(time_PR < timestamp(1)-.5 & time_PR > session_PPG_time(1)+.5);
 
    if id > 200
        HRV_supine_old = calculate_HRVfeatures(RR_win*60,time_RR_win*60,Fs_resample,HRV_supine_old);
        PRV_supine_old = calculate_HRVfeatures(PR_win*60,time_PR_win*60,Fs_resample,PRV_supine_old);
    else
        HRV_supine_young = calculate_HRVfeatures(RR_win*60,time_RR_win*60,Fs_resample,HRV_supine_young);
        PRV_supine_young = calculate_HRVfeatures(PR_win*60,time_PR_win*60,Fs_resample,PRV_supine_young);       
    end
    %% upright
    
    RR_win = RR(time_RR < session_PPG_time(end)-.5 & time_RR > timestamp(2)+.5);
    PR_win = PR(time_PR < session_PPG_time(end)-.5 & time_PR > timestamp(2)+.5);
    time_RR_win = time_RR(time_RR < session_PPG_time(end)-.5 & time_RR > timestamp(2)+.5);
    time_PR_win = time_PR(time_PR < session_PPG_time(end)-.5 & time_PR > timestamp(2)+.5);

    if id > 200
        HRV_upright_old = calculate_HRVfeatures(RR_win*60,time_RR_win*60,Fs_resample,HRV_upright_old);
        PRV_upright_old = calculate_HRVfeatures(PR_win*60,time_PR_win*60,Fs_resample,PRV_upright_old);    
    else
        HRV_upright_young = calculate_HRVfeatures(RR_win*60,time_RR_win*60,Fs_resample,HRV_upright_young);
        PRV_upright_young = calculate_HRVfeatures(PR_win*60,time_PR_win*60,Fs_resample,PRV_upright_young);
    end    
    %% transition
    
    RR_win = RR(time_RR < timestamp(2)+1 & time_RR > timestamp(1)-1);
    PR_win = PR(time_PR < timestamp(2)+1 & time_PR > timestamp(1)-1);
    time_RR_win = time_RR(time_RR < timestamp(2)+1 & time_RR > timestamp(1)-1);
    time_PR_win = time_PR(time_PR < timestamp(2)+1 & time_PR > timestamp(1)-1);

    if id > 200
        HRV_transition_old = calculate_HRVfeatures(RR_win*60,time_RR_win*60,Fs_resample,HRV_transition_old);
        PRV_transition_old = calculate_HRVfeatures(PR_win*60,time_PR_win*60,Fs_resample,PRV_transition_old);
    else
        HRV_transition_young = calculate_HRVfeatures(RR_win*60,time_RR_win*60,Fs_resample,HRV_transition_young);
        PRV_transition_young = calculate_HRVfeatures(PR_win*60,time_PR_win*60,Fs_resample,PRV_transition_young);
    end
end

%% statistical analysis HRV-PRV

% Young 
results_transitions_young = consistent_analysis(HRV_transition_young,PRV_transition_young);
results_supine_young = consistent_analysis(HRV_supine_young,PRV_supine_young);
results_upright_young = consistent_analysis(HRV_upright_young,PRV_upright_young);
% 
% % Old
% results_transitions_old = consistent_analysis(HRV_transition_old,PRV_transition_old);
% results_supine_old = consistent_analysis(HRV_supine_old,PRV_supine_old);
% results_upright_old = consistent_analysis(HRV_upright_old,PRV_upright_old);


%% plots 
data_HRV.HRV_transition_old = HRV_transition_old; data_HRV.HRV_transition_young = HRV_transition_young;
data_HRV.HRV_supine_old = HRV_supine_old; data_HRV.HRV_supine_young = HRV_supine_young;
data_HRV.HRV_upright_old = HRV_upright_old;  data_HRV.HRV_upright_young = HRV_upright_young;
data_PRV.PRV_transition_old = PRV_transition_old; data_PRV.PRV_transition_young = PRV_transition_young;
data_PRV.PRV_supine_old = PRV_supine_old; data_PRV.PRV_supine_young = PRV_supine_young;
data_PRV.PRV_upright_old = PRV_upright_old;  data_PRV.PRV_upright_young = PRV_upright_young;

plot_figures_hrvprv(data_HRV,data_PRV,'young') % HRV-PRV

