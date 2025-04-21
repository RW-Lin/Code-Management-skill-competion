function output = calculate_HRVfeatures(RR_interval,time_idx,Fs_resample,HRV_all)
%% 
    % RR_interval: in seconds
    % RR_interval = RR_interval*1e3;
    output = HRV_all;
    
%% Time domain
    output.SDNN = [output.SDNN 1e3*HRV.SDNN(RR_interval,0,1)];
    output.RMSSD = [output.RMSSD 1e3*HRV.RMSSD(RR_interval,0,1)];
    output.pNN50 = [output.pNN50 100*HRV.pNN50(RR_interval,0,1)];
    output.AHR = [output.AHR HRV.HR(RR_interval,0)];    
    output.SDSD = [output.SDSD 1e3*HRV.SDSD(RR_interval,0,1)]; 
    % output.TINN = [output.TINN HRV.TINN(RR_interval,0)];
    % output.TRI = [output.TRI HRV.TRI(RR_interval)];
%% Frequency domain
    % [pLF, pHF, LFHFratio, VLF, LF, HF] = HRV.fft_val_fun(RR_interval, Fs_resample, 'pchip'); 
    [LF,HF,pLF,pHF,LFHFratio]  = calc_hrv_comparison(time_idx,RR_interval*1e3);
    output.pLF = [output.pLF pLF]; output.pHF = [output.pHF pHF];
    output.LFHFratio = [output.LFHFratio LFHFratio]; 
    output.LF = [output.LF LF]; output.HF = [output.HF HF];
    % output.TotalP = [output.LF TotalPower];
    
%% Nonlinear domain

    [SD1, SD2, SD1SD2ratio] = HRV.returnmap_val(RR_interval, 0);
    % output.ApEn = [output.ApEn HRV.ApEn(RR_interval,0)];
    output.SD1 = [output.SD1 1e3*SD1];  output.SD2 = [output.SD2 1e3*SD2];
    output.SDratio = [output.SDratio SD1SD2ratio];     

    % entropy
    % interp RR series
    % valid_idx =  ~isnan(RR_interval(:)); 
    % new_time_RR = min(time_idx(valid_idx)):1/Fs_resample:max(time_idx(valid_idx));
    % RR_intp = interp1(time_idx(valid_idx),RR(valid_idx)*60,new_time_RR,'pchip');
    
    lag = 1; dim = 2;
    ApEn = approximateEntropy(RR_interval(~isnan(RR_interval(:))), lag, dim); %Approximate entropy with default 
    SampEn = sampen(RR_interval(~isnan(RR_interval(:))), dim, .2); % sample entropy
    output.ApEn = [output.ApEn ApEn];  output.SampEn = [output.SampEn SampEn];
end

function [LF_power,HF_power,LFP_fft,HFP_fft,LF_HF_ratio_fft] =...
    calc_hrv_comparison(time_data, rr_data)
% calc_hrv_comparison 对比 Lomb-Scargle 和插值+PSD 的 HRV 频谱计算

%% Step 1: 数据清理（去除 NaN）
valid_idx = ~isnan(rr_data);
time_clean = time_data(valid_idx);
rr_clean = rr_data(valid_idx);


%% Step 3: 插值 + PSD 计算
target_fs = 4;  % 4Hz 采样率
dt = 1 / target_fs;
time_interp = min(time_clean):dt:max(time_clean);

rr_interp = interp1(time_clean, rr_clean, time_interp, 'pchip');
N = length(rr_interp);  % 信号长度
Y = fft(rr_interp);  % 快速傅里叶变换
P2 = abs(Y / N).^2;  % 双边功率谱密度
P1 = P2(1:floor(N/2)+1);  % 单边功率谱
P1(2:end-1) = 2 * P1(2:end-1);  % 单边谱补偿
f_fft = (0:floor(N/2)) * (target_fs / N);  % 频率向量

% 计算低频 (LF) 和 高频 (HF) 功率
LF_band = [0.04, 0.15];  % LF 频率范围 (Hz)
HF_band = [0.15, 0.4];   % HF 频率范围 (Hz)

LF_power = sum(P1(f_fft >= LF_band(1) & f_fft <= LF_band(2)));  % LF 功率
HF_power = sum(P1(f_fft >= HF_band(1) & f_fft <= HF_band(2)));  % HF 功率
total_power = sum(P1(f_fft >= LF_band(1) & f_fft <= HF_band(2)));  % 总功率 (LF + HF)

% 计算 LFP, HFP 和 LF/HF 比值
LFP_fft = LF_power / total_power;
HFP_fft = HF_power / total_power;
LF_HF_ratio_fft = LF_power / HF_power;

end

