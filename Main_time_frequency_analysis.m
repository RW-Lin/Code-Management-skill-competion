% Main time frequency analysis HRV/PRV/PAT
% compute coherence (t), spearman CC of HF & LF & LHFHF ratio
close all, clear all,
Fs_resample = 4; %s^-1
lag_before = 5; lag_after = 3; % duration before and after timestamp
fs_ppg = 1e3; debug = 0;
activity = {'Supine_stand','HUT','SittoStand','VM','RapidSupinestand','SlowSupinestand'};  
[b,a] = butter(5,[0.01,0.8]/(Fs_resample/2));

%% test
for act = 1:length(activity)
savepath = ['C:\Users\LinR\OneDrive - University of Twente\work\project\Spectral analysis of HRV, PTT and BP\Codes\store_results\' activity{act} '\' ]; 
dirInfo = dir(savepath);
LFcc_all = zeros(length(dirInfo)-2,1); HFcc_all = zeros(length(dirInfo)-2,1); LFHFcc_all = zeros(length(dirInfo)-2,1);
% rhos_all = zeros(length(dirInfo)-2,1);
for sg_num = 3 %3:length(dirInfo) % 27 37 45 48

load([savepath dirInfo(sg_num).name])
data_ppg = data_segments.ppg; data_ecg = data_segments.ecg;
ppg_peaktime = data_segments.ppgpeak_time; ecg_peaktime = data_segments.ecgpeak_time;
session_PPG_time = data_segments.ppgtime; session_ECG_time = data_segments.ecgtime;
ppg_len = length(data_segments.ppg) / 1e3; timestamp = data_segments.timestamp;
psqi_session = data_segments.psqi; 
ecg_peakidx = data_segments.ecg_idx; ppg_peakidx = data_segments.ppg_idx;
[PAT, calibrated_peaks_ecg, calibrated_peaks_ppg] =...
        calculate_PATv(ecg_peaktime*60,ppg_peaktime*60,psqi_session,fs_ppg,timestamp*60,lag_before);

RR = HRV.RRfilter(diff(calibrated_peaks_ecg),60); 
PR = HRV.RRfilter(diff(calibrated_peaks_ppg),60);

time_RR = calibrated_peaks_ecg(1:end-1)*60; time_PR = calibrated_peaks_ppg(1:end-1)*60;
PAT(PAT<100 | PAT>500) = NaN; time_PAT = calibrated_peaks_ppg;

% HRV
valid_idx =  ~isnan(RR(:)) & ~isnan(PR(:)) & ~isnan(time_RR(:)) & ~isnan(time_PR(:)); 
new_time_RR = min(time_RR(valid_idx)):1/Fs_resample:max(time_RR(valid_idx));
RR_intp = interp1(time_RR(valid_idx),RR(valid_idx)*60,new_time_RR,'pchip');
% RR_intp = filtfilt(b,a,RR_intp);
% PRV
new_time_PR = min(time_PR(valid_idx)):1/Fs_resample:max(time_PR(valid_idx));
PR_intp = interp1(time_PR(valid_idx),60*PR(valid_idx),new_time_PR,'pchip');
% PR_intp = filtfilt(b,a,PR_intp);
% valid_idx = ~isnan(HRV.RRfilter(PR_intp(:),60)) & ~isnan(HRV.RRfilter(RR_intp(:),60)); 
% RR_intp = RR_intp(valid_idx); PR_intp = PR_intp(valid_idx);

if debug
figure,
subplot(411)
plot(session_ECG_time,data_ecg);hold on;plot(ecg_peaktime,data_ecg(ecg_peakidx),'o')
plot(session_PPG_time,data_ppg);hold on;plot(ppg_peaktime,data_ppg(ppg_peakidx),'o')
xline(timestamp)
subplot(412)
plot(time_RR,RR,'*-');hold on;plot(new_time_RR,RR_intp,'s-')
plot(time_PR,PR,'*-');hold on;plot(new_time_PR,PR_intp,'s-')
xline(timestamp*60)
subplot(413)
imagesc(t, f, (abs(tfr_hrv)));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('HRV Time-Frequency Representation (Pseudo Wigner-Ville)');
colorbar;
subplot(414)
imagesc(t, f, abs(tfr_prv));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('PRV Time-Frequency Representation (Pseudo Wigner-Ville)');
colorbar;
end

% remove possible extra RR-/PR interval
if length(RR_intp)-length(PR_intp) == 1
    RR_intp(end) = [];
elseif length(RR_intp) - length(PR_intp) == -1
    PR_intp(end) = [];
end
[Metrics, HRV_TF, PRV_TF,t,f] = time_frequency_analysis(RR_intp,PR_intp,Fs_resample);
LFHFcc_all(sg_num-2) = Metrics.LFHFcc; LFcc_all(sg_num-2) = Metrics.LFcc; HFcc_all(sg_num-2) = Metrics.HFcc;
% rhos_all = Metrics.rho_s;
end
figure,boxplot([LFHFcc_all,LFcc_all,HFcc_all],'Labels',{'Ratio','LF','HF'})
end









%%

%%
% try spwvd
% valid_idx = ~isnan(PAT(:)) & ~isnan(time_PAT(:)); new_time_PAT = min(time_PAT)*60:1/Fs_resample:60*max(time_PAT);
% PAT_intp = interp1(time_PAT(valid_idx)*60,PAT(valid_idx),new_time_PAT,'pchip');
% PAT_intp = filtfilt(b,a,PAT_intp);
% % [tfr,f,t] = wvd(RR_intp,Fs_resample,"smoothedPseudo",kaiser(261),kaiser(121),'NumFrequencyPoints',ceil(length(RR_intp)/2));
% [tfr_pat,f,t] = wvd(PAT_intp,Fs_resample,"smoothedPseudo");
% figure,
% subplot(311)
% plot(time_PAT,PAT);xline(timestamp)
% subplot(312)
% plot(time_PAT*60,PAT);hold on;plot(new_time_PAT,PAT_intp,'s-')
% xline(timestamp*60)
% subplot(313)
% imagesc(t, f, log10(abs(tfr_pat)));
% axis xy;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('HRV Time-Frequency Representation (Pseudo Wigner-Ville)');
% colorbar;
