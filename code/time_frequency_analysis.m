function [Metrics, HRV_TF, PRV_TF,t,f] = time_frequency_analysis(hrv_signal,prv_signal,Fs_resample)
%% initialization

%% interploration


%% perform MSTP
% [t, f, HRV_TF] = multitaper_spectrogram(hrv_signal, fs, time_window);
% [~, ~, PRV_TF] = multitaper_spectrogram(prv_signal, fs, time_window);
%%
% t = 0:1/Fs_resample:(length(hrv_signal)-1)/Fs_resample;
% % N = 2^nextpow2(length(hrv_signal));
% % glength = floor(N/20)+1; hlength = floor(N/20)+1;
% % 
% % g = gausswin(glength); h = gausswin(hlength);
% % 
% % [tfxy,~,~] = tfrspwv([hrv_signal(:),prv_signal(:)],1:length(hrv_signal),N,g,h);
% % [tf_hrv,t,f] = tfrspwv(hrv_signal(:),1:length(hrv_signal),N,g,h);
% % [tf_prv,~,f] = tfrspwv(prv_signal(:),1:length(hrv_signal),N,g,h);
% % 
% % [N_freq, ~] = size(tfxy); % 获取时频分布的大小
% % tfr_hrv_positive = tf_hrv(1:N_freq/2, :); % 提取正频率部分
% % tfr_prv_positive = tf_prv(1:N_freq/2, :); % 提取正频率部分
% % 
% % % % 频率向量
% % N_freq = size(tfxy, 1); f = f*Fs_resample;
% % tfr_positive = tfxy(1:N_freq/2, :) * 2;  % 正频率部分乘以 2
% % f_positive = f(1:N_freq/2);
% % 
% % % 绘图
% % figure,
% % subplot(311)
% % imagesc(t/Fs_resample, f_positive, abs(tfr_positive));
% % xlabel('Time (s)');
% % ylabel('Frequency (Hz)');
% % title('Positive Frequency SPWVD');
% % axis xy; colorbar;
% % subplot(312)
% % imagesc(t/Fs_resample, f_positive, abs(tfr_hrv_positive));
% % xlabel('Time (s)');
% % ylabel('Frequency (Hz)');
% % title('Positive Frequency SPWVD');
% % axis xy; colorbar;
% % subplot(313)
% % imagesc(t/Fs_resample, f_positive, abs(tfr_prv_positive));
% % xlabel('Time (s)');
% % ylabel('Frequency (Hz)');
% % title('Positive Frequency SPWVD');
% % axis xy; colorbar;
% % gamma2_band = coherence_band(tfr_hrv_positive, tfr_prv_positive, tfr_positive, f_positive, band)
%% SPWVD
% [HRV_TF, t, f] = multitaper_spectrogram(hrv_signal(:), Fs_resample);
% [PRV_TF, t, f] = multitaper_spectrogram(prv_signal(:), Fs_resample);

twin_length = round(400 / 1.25)+1; fwin_length = round(400 / 1.25)+1;   
twin = gausswin(twin_length); fwin = gausswin(fwin_length);
% 
% [tfxy,t2,f2] = tfrspwv([hrv_signal(:),prv_signal(:)],1:length(hrv_signal));
[TFxy,f1,t1] = xwvd(hrv_signal,prv_signal,Fs_resample,"smoothedPseudo");
[HRV_TF,f,t] = wvd(hrv_signal,Fs_resample,"smoothedPseudo");
[PRV_TF,~,~] = wvd(prv_signal,Fs_resample,"smoothedPseudo");
% figure;
% subplot(4,1,1);
% imagesc(t,f,(HRV_TF)); axis xy; title('HRV TF Magnitude');colorbar;axis xy
% subplot(4,1,2);
% imagesc(t,f,(PRV_TF)); axis xy; title('PRV TF Magnitude');colorbar;axis xy
% subplot(4,1,3);
% imagesc(t1,f1,(TFxy)); axis xy; title('Cross TF Magnitude');colorbar;axis xy


%% Metrics
LF_band = [0.04 0.15]; HF_band = [0.15 0.4];

HRV_LF = integrate_band(HRV_TF, f, LF_band); PRV_LF = integrate_band(PRV_TF, f, LF_band);
LFcc = corr(PRV_LF(:),HRV_LF(:),'type','Spearman');
HRV_HF = integrate_band(HRV_TF, f, HF_band); PRV_HF = integrate_band(PRV_TF, f, HF_band);
HFcc = corr(PRV_HF(:),HRV_HF(:),'type','Spearman');
HRV_LFHFratio = HRV_LF./HRV_HF; PRV_LFHFratio = PRV_LF./PRV_HF;
LFHFcc = corr(PRV_LFHFratio(:),HRV_LFHFratio(:),'type','Spearman');

% energy
delta_LF = HRV_LF - PRV_LF;
delta_HF = HRV_HF - PRV_HF;
% 
figure,
subplot(511)
plot(60./hrv_signal);hold on;plot(60./prv_signal);legend('HRV','PRV')
subplot(512)
plot(HRV_HF);hold on;plot(PRV_HF);ylabel('HF')
subplot(513)
plot(HRV_LF);hold on;plot(PRV_LF);
legend('HRV-LF','PRV-LF'); ylabel('LF')
subplot(514)
plot(HRV_LFHFratio); hold on;plot(PRV_LFHFratio);ylabel('LF-HF ratio')
legend('HRV','PRV')
subplot(515)
plot(delta_HF); hold on; plot(delta_LF);legend('HF','LF')
% ******************************************
% band coherence
time_window_size = 4; % 平滑窗口（8 个时间点）
[Coherence_t_f, t_coh, f_coh] = compute_coherence_from_tfr(PRV_TF, HRV_TF, TFxy, f, t, time_window_size);

LF_idx = (f_coh >= LF_band(1)) & (f_coh <= LF_band(2));
HF_idx = (f_coh >= HF_band(1)) & (f_coh <= HF_band(2));

% 平均频率带内的相干性随时间变化
Coherence_LF_t = mean(Coherence_t_f(LF_idx, :), 1); % LF 带
Coherence_HF_t = mean(Coherence_t_f(HF_idx, :), 1); % HF 带

% 可视化频带相干性随时间变化
figure;
plot(t_coh, Coherence_LF_t, 'LineWidth', 2);
hold on;
plot(t_coh, Coherence_HF_t, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Coherence');
title('LF and HF Band Coherence over Time');
legend('LF Band', 'HF Band');
grid on;
% correlation
% rho_s = tf_corr(HRV_TF, PRV_TF, t);



% Coherence
% gamma = tf_coherence(HRV_TF, PRV_TF);


% 
% Metrics.rho_s = rho_s;
Metrics.PRV_LF=PRV_LF; Metrics.PRV_LF=PRV_LF;
Metrics.HRV_LF=HRV_LF; Metrics.HRV_HF=HRV_HF;
% Metrics.delta_LF = delta_LF; Metrics.delta_HF = delta_HF;
% Metrics.gamma_LF = gamma2_LF; Metrics.gamma_HF = gamma2_HF;
Metrics.LFcc = LFcc; Metrics.HFcc = HFcc;
Metrics.LFHFcc = LFHFcc; 
%%
debug=0;
if debug
    figure,
    subplot(311)
    plot(HRV_LF);hold on;plot(PRV_LF);title('HF')
    subplot(312)
    plot(HRV_HF);hold on;plot(PRV_HF); legend('HRV','PRV');legend boxoff; title('HF')
    subplot(313)
    plot(delta_LF);hold on;plot(delta_HF);title('Energy difference');legend('LF','HF')

    % figure,
    subplot(511)
    subplot(512)
    subplot(513)
    plot(t,rho_s);title('spectral correlation')
    subplot(514)
    plot(t,delta_LF);hold on;plot(t,delta_HF);legend('\DeltaLF','\DeltaHF')
    subplot(515)
    imagesc(t,f,gamma);
    colorbar
end
end

function power_band = integrate_band(TF, f, band)
    idx = (f >= band(1)) & (f <= band(2));  % 选择频率索引
    power_band = sum(abs(TF(idx, :).^2), 1); % 功率密度积分
end

function gamma2 = tf_coherence(TF1, TF2)
% Time-frequency coherence calculation (pointwise)

% Ensure TF1 and TF2 have the same size
if ~isequal(size(TF1), size(TF2))
    error('TF1 and TF2 must have the same dimensions.');
end

% Numerator: |TF1 * conj(TF2)|^2 (element-wise)
num = abs(TF1 .* conj(TF2)).^2;

% Denominator: |TF1|^2 * |TF2|^2 (element-wise)
den = (abs(TF1).^2) .* (abs(TF2).^2);

% Coherence calculation (element-wise division)
gamma2 = num ./ den;

end
function [Coherence_t_f, t_coh, f_coh] = compute_coherence_from_tfr(TFR_PRV, TFR_HRV, TFxy, f, t, time_window_size)
% compute_coherence_from_tfr - 基于已计算的 TFR 计算 HRV 和 PRV 的相干性
%
% 输入参数：
% TFR_PRV          - PRV 的时频分布 (矩阵大小为 [频率点, 时间点])
% TFR_HRV          - HRV 的时频分布 (矩阵大小为 [频率点, 时间点])
% f                - 频率向量 (Hz)
% t                - 时间向量 (秒)
% time_window_size - 相干性计算的时间平滑窗口大小（单位：时间点）
%
% 输出参数：
% Coherence_t_f    - 相干性随时间和频率的二维矩阵
% t_coh            - 时间向量
% f_coh            - 频率向量

    % 1. 检查输入大小
    assert(size(TFR_PRV, 1) == size(TFR_HRV, 1), '频率维度不匹配');
    assert(size(TFR_PRV, 2) == size(TFR_HRV, 2), '时间维度不匹配');

    % 2. 自谱密度和交叉谱密度
    AutoSpectrum_PRV = abs(TFR_PRV).^2; % PRV 自谱
    AutoSpectrum_HRV = abs(TFR_HRV).^2; % HRV 自谱
    CrossSpectrum = abs(TFxy).^2; % 交叉谱

    % 3. 平滑自谱和交叉谱密度
    smooth_func = @(x) movmean(x, time_window_size, 2); % 时间维度平滑
    Smooth_Auto_PRV = smooth_func(AutoSpectrum_PRV);
    Smooth_Auto_HRV = smooth_func(AutoSpectrum_HRV);
    Smooth_Cross = smooth_func(CrossSpectrum);

    % 4. 计算相干性
    Coherence_t_f = Smooth_Cross ./ (Smooth_Auto_PRV .* Smooth_Auto_HRV);

    % 5. 时间和频率向量
    t_coh = t; % 时间向量直接沿用输入
    f_coh = f; % 频率向量直接沿用输入
end


function rho_s = spectral_correlation(TF1, TF2)
% Instantaneous spectral correlation
rho_s = zeros(1, size(TF1, 2));
for t = 1:size(TF1, 2)
    rho_s(t) = corr(TF1(:, t), TF2(:, t),'type','Spearman');
end
end

function rho_s = tf_corr(TF1, TF2, t)
% TF_CORR: 计算两个信号在给定时间点 t0 的时频谱的 Pearson 相关系数
%
% 输入:
%   TF1 : 第一个信号的时频表示 (例如 WVD)
%   TF2 : 第二个信号的时频表示 (同尺寸)
%   t   : 时间向量
%   t0  : 指定的时间点
%
% 输出:
%   rho_s : 在时间点 t0 两个信号时频谱的 Pearson 相关系数

    % 确保 TF1 和 TF2 维度一致
    if ~isequal(size(TF1), size(TF2))
        error('TF1 和 TF2 的维度必须一致');
    end
    
    rho_s = zeros(size(t));
    for t_idx = 1:length(t)
        % 提取时间点 t0 的频谱
        S1 = TF1(:, t_idx);  % 第一个信号的时频谱
        S2 = TF2(:, t_idx);  % 第二个信号的时频谱
        
        % 计算 Pearson 相关系数
        rho_s(t_idx) = corr(S1, S2, 'Type', 'Spearman');
    end
end