function [output,id,count] = processing_posturalchange(data,annotation,name)
close all
savepath = 'C:\Users\LinR\OneDrive - University of Twente\work\project\Spectral analysis of HRV, PTT and BP\Codes\store_results\HGS\';
% This function only evaluates Supine2Upright
    %% initialization
    fs_ppg = 1e3; Fs_resample = 4; fields = fieldnames(annotation);
    [b_ppg,a_ppg] = butter(2,[.67 8]/(fs_ppg/2));
    [b,a] = butter(4,0.003/(Fs_resample/2),'high');
    debug=0; % control plots
    lag_before=5; lag_after = 3; % minutes before/after test 

    good_list = [];
    count = 0;
    options_ecg.beat_detectors = {'jqrs'}; beat_detector = 'msptdfastv2';
    if eval(name)>199, id = 0; else, id = 1; end % id (annann0->old 1->young) 
    frqrange = [0.01 0.8]; taperparams = [2 3]; windowparams = [26.5,5]; % parameters for MTSP
    
    win_len = 120; % second 
    overlap = win_len * 0.5; 
    %% BP
    SBP_time = data.reSYS.time; data_SBP = data.reSYS.values; 

    %% ECG
    fs_ecg = 300;  ECG = data.ECG_aVF.values; ECG_time = data.ECG_aVF.time; 
    [b_ecg,a_ecg] = butter(3,[.5 120]/(fs_ecg/2)); 
    ecg_data = filtfilt(b_ecg,a_ecg,ECG); 
    % ecg_data = ECG;

    %% PPG
    PPG_finger = data.PPG2.values; PPG_time = data.PPG2.time; 
    PPG_finger(isnan(PPG_finger)) = 0;
    PPG_finger = -filtfilt(b_ppg,a_ppg,PPG_finger);
    s.v=PPG_finger; s.fs=fs_ppg; 
    
    %%  timestamp
    if id == 0,
        if isfield(annotation,"HGS") 
           time_stamps = [annotation.HGS];
        % elseif isfield(annotation,"Rapid_supine_to_stand")
        %    time_stamps = [annotation.Rapid_supine_to_stand;annotation.Slow_supine_to_stand];
        end
    % else
    %     if isfield(annotation,"Supine_to_stand")
    %         time_stamps = [annotation.Supine_to_stand];
    %     elseif isfield(annotation,"Rapid_supine_to_stand")
    %        time_stamps = [annotation.Rapid_supine_to_stand;annotation.Slow_supine_to_stand];
    %     end
    end
    remove_id = [];
    
    for ts = 1:size(time_stamps,1) 
        if time_stamps(ts,1) > max(PPG_time) || time_stamps(ts,1) < min(PPG_time)
            remove_id = [remove_id ts];
        elseif time_stamps(ts,2) > max(PPG_time) || time_stamps(ts,2) <  min(PPG_time)
            remove_id = [remove_id ts];
        elseif time_stamps(ts,1) > ECG_time(end) || time_stamps(ts,1) <  min(PPG_time)
            remove_id = [remove_id ts];
        elseif time_stamps(ts,2) > ECG_time(end) || time_stamps(ts,2) <  min(PPG_time)
            remove_id = [remove_id ts];
        end
    end
    time_stamps(unique(remove_id),:) = [];
   
    %% analyze session data
    % In this for loop we compute coherence for each session
    
    % HRV_all.Supine.RMSSD = zeros(size(time_stamps,1),(lag_before*60-(segment_len))/10+1);
    % HRV_all.Upright.SDNN = zeros(size(time_stamps,1),(lag_after*60-(segment_len))/10+1);
    % HRV_all.Supine.pNN50 = zeros(size(time_stamps,1),(lag_after*60-(segment_len))/10+1);
    % HRV_all.Upright.pNN50 = zeros(size(time_stamps,1),(lag_after*60-(segment_len))/10+1);    
    % PRV_all.Supine.RMSSD = zeros(size(time_stamps,1),(lag_before*60-(segment_len))/10+1);
    % PRV_all.Upright.SDNN = zeros(size(time_stamps,1),(lag_after*60-(segment_len))/10+1);

    for session = 1:size(time_stamps,1)
        
        timestamp = time_stamps(session,:);

        % session time & data
        session_index_ppg = find(PPG_time>timestamp(1)-lag_before & PPG_time<timestamp(2)+lag_after);
        data_ppg_session = PPG_finger(session_index_ppg); session_PPG_time = PPG_time(session_index_ppg);

        session_index_ecg = find(ECG_time>=timestamp(1)-lag_before & ECG_time<=timestamp(2)+lag_after);
        data_ecg_session = ecg_data(session_index_ecg); session_ECG_time = ECG_time(session_index_ecg);

        session_index_sbp = find(SBP_time>=timestamp(1)-lag_before & SBP_time<=timestamp(2)+lag_after);
        %% 
        data_sbp_session = data_SBP(session_index_sbp); session_SBP_time = SBP_time(session_index_sbp);
        
        % detect peaks
        peaks_ecg_session = detect_ecg_peaks(data_ecg_session,fs_ecg);  
        
        s.v=data_ppg_session(:); s.fs=fs_ppg; 
        [~, ~, peaks_ppg_session] = detect_ppg_beats(s, beat_detector);
        peaks_ecg_time = session_ECG_time(peaks_ecg_session); peaks_ppg_time = session_PPG_time(peaks_ppg_session);
        
        % R-R interval and PPG signal quality index
        % RR_interval = diff(peaks_ecg_time*60); PR_interval = diff(peaks_ppg_time*60);
        % HR = 60./RR_interval; PR = 60./PR_interval;
        psqi_session=[ppgSQI(data_ppg_session,peaks_ppg_session,fs_ppg) 0.5];


        % select data for analysis
        [accept_idx,gd_ratio] = check_ppg_sqi(psqi_session(1:end-1));
        if  accept_idx && abs(length(peaks_ecg_time)-length(peaks_ppg_time))<20 ...
            && sum(peaks_ppg_time<timestamp(1))>200 ...
            && sum(peaks_ppg_time>timestamp(2))>100 && sum(peaks_ecg_time<timestamp(1))>200 ...
            && sum(peaks_ecg_time>timestamp(2))>100 && mean(data_ecg_session(peaks_ecg_session))>1e-2
            
            good_list = [good_list session]; 
   
            % figure, subplot(311),
            % plot(session_PPG_time*60,data_ppg_session); ylabel('PPG [n.u.]')
            % hold on; plot(peaks_ppg_time*60,data_ppg_session(peaks_ppg_session),'o')
            % xline([timestamp(1)*60 timestamp(2)*60],'g','LineWidth',1.5);
            % subplot(312); plot(session_ECG_time*60,data_ecg_session);
            % xline([timestamp(1)*60 timestamp(2)*60],'g','LineWidth',1.5);
            % hold on; plot(peaks_ecg_time*60,data_ecg_session(peaks_ecg_session),'o');
            % xline([timestamp(1)*60 timestamp(2)*60],'g','LineWidth',1.5);ylabel('ECG [mV]')
            % subplot(313)
            % plot(peaks_ecg_time(1:end-1)*60,HR,'*-','LineWidth',2);hold on;
            % plot(peaks_ppg_time(1:end-1)*60,PR,'*-','LineWidth',1.25);
            % xline([timestamp(1)*60 timestamp(2)*60],'g','LineWidth',2);
            % legend('HR','PR'); ylabel('[bpm]');xlabel('Time [s]')
            % PAT
            % [PAT, calibrated_peaks_ppg, calibrated_peaks_ecg] =...
            %     calculate_PATv(peaks_ecg_time*60,peaks_ppg_time*60,psqi_session,fs_ppg,timestamp*60,lag_before);
            data_segments.psqi = psqi_session;data_segments.fsppg = fs_ppg; data_segments.ratio = gd_ratio;  
            data_segments.fs_ecg = fs_ecg; data_segments.timestamp = timestamp;
            data_segments.ppg = data_ppg_session; data_segments.ecg = data_ecg_session;
            data_segments.ecg_idx = peaks_ecg_session; data_segments.ppg_idx = peaks_ppg_session;
            data_segments.ppgpeak_time = peaks_ppg_time; data_segments.ecgpeak_time = peaks_ecg_time;
            data_segments.ppgtime = session_PPG_time; data_segments.ecgtime = session_ECG_time;
            data_segments.sbptime = session_SBP_time; data_segments.sbp = data_sbp_session;
            save([savepath 'Data_' name '_session_' num2str(session) ],"data_segments")
            % 
            % HRV_features = calculate_HRVfeatures(calibrated_peaks_ecg,Fs_resample,HRV_features_all); 
            % PRV_features = calculate_HRVfeatures(calibrated_peaks_ppg,Fs_resample,PRV_features_all);                    
        end
        
        % PAT
        
    %%
        

        % plot HR vs. PR
        % figure, 
        % subplot(411),plot(1/fs_ppg:1/fs_ppg:length(data_ppg_session)/fs_ppg,data_ppg_session);
        % hold on; plot(peaks_ppg_session/fs_ppg,data_ppg_session(peaks_ppg_session),'o')
        % plot(peaks_bp_session/fs_ppg,data_ppg_session(peaks_bp_session),'*')
        % xline([60*lag_before 60*lag_before+diff(timestamp)*60],'g','LineWidth',2);ylabel("PPG")
        % subplot(412),
        % plot(1/fs_ppg:1/fs_ppg:length(data_ecg_session)/fs_ppg,data_ecg_session);
        % hold on; plot(peaks_ecg_session/fs_ppg,data_ecg_session(peaks_ecg_session),'o');
        % xline([60*lag_before 60*lag_before+diff(timestamp)*60],'g','LineWidth',2);ylabel("ECG")
        % subplot(413)
        % plot(peaks_ecg_session(1:end-1)/fs_ppg,HR,'*-','LineWidth',2);hold on;
        % plot(peaks_ppg_session(1:end-1)/fs_ppg,PR,'*-');
        % % plot(peaks_bp_session(1:end-1)/fs_ppg,BPR,'*-');
        % ylabel("Raw Rate [bpm]"); xline([60*lag_before 60*lag_before+diff(timestamp)*60],'g','LineWidth',1.5)
        % legend('HR','PR','tilt start','tilt end');
        % subplot(414), plot(new_time_stamphr,HR_resample,'*-','LineWidth',2);hold on;
        % plot(new_time_stamppr,PR_resample,'*-');ylabel("Resampled Rate [bpm]")
        % xline([60*lag_before 60*lag_before+diff(timestamp)*60],'g','LineWidth',1.5)
        % sgtitle(['Session' num2str(session) ' - Subject-' name])
        
    % end

    
    %% before test / upright
    % for sg = 1:(lag_before*60-(segment_len))/overlap+1
    %     time_start = 1+(sg-1)*overlap*fs_ppg; time_end = (sg-1)*overlap*fs_ppg+segment_len*fs_ppg;
    %     RR_seg = RR_interval(peaks_ecg_session(1:end-1)>time_start & peaks_ecg_session(1:end-1)<time_end);
    %     PR_seg = PR_interval(peaks_ppg_session(1:end-1)>time_start & peaks_ppg_session(1:end-1)<time_end);
    %     psqi_sg = psqi_session(peaks_ppg_session(1:end-1)>time_start & peaks_ppg_session(1:end-1)<time_end);
    %     if length(RR_seg)-length(PR_seg)>10 || check_ppg_sqi(psqi_sg)==0,  continue; end 
    %     % figure,plot(RR_seg,'*-');hold on;plot(PR_seg,'o-'); legend('HR','PR')
    %     % compute features
    %     SDNN_HR = HRV.SDNN(RR_seg,0,1); SDNN_PR = HRV.SDNN(PR_seg,0,1);
    %     RMSSD_HR = HRV.RMSSD(RR_seg,0,1); RMSSD_PR = HRV.RMSSD(PR_seg,0,1);
    %     pNN50_HR = HRV.pNN50(RR_seg,0,1); pNN50_PR = HRV.pNN50(PR_seg,0,1);
    %     ApEn_HR = HRV.ApEn(RR_seg); ApEn_PR = HRV.ApEn(PR_seg);
    % 
    %     HRV_all.Supine.RMSSD(session,sg) = RMSSD_HR; HRV_all.Supine.SDNN(session,sg) = SDNN_HR;
    %     PRV_all.Supine.RMSSD(session,sg) = RMSSD_PR; PRV_all.Supine.SDNN(session,sg) = SDNN_PR;
    %     HRV_all.Supine.pNN50(session,sg) = pNN50_HR; PRV_all.Supine.pNN50(session,sg) = pNN50_PR;
    %     HRV_all.Supine.ApEn(session,sg) = ApEn_HR; PRV_all.Supine.ApEn(session,sg) = ApEn_PR;
    % end
    % % after / upright
    % for sg = 1:(lag_after*60-(segment_len))/overlap+1
    %     time_start = 1+(sg-1)*overlap*fs_ppg+(60*lag_before+diff(timestamp))*fs_ppg;
    %     time_end = (sg-1)*overlap*fs_ppg+segment_len*fs_ppg+(60*lag_before+diff(timestamp))*fs_ppg;
    %     RR_seg = RR_interval(peaks_ecg_session(1:end-1)>time_start & peaks_ecg_session(1:end-1)<time_end);
    %     PR_seg = PR_interval(peaks_ppg_session(1:end-1)>time_start & peaks_ppg_session(1:end-1)<time_end);
    %     psqi_sg = psqi_session(peaks_ppg_session(1:end-1)>time_start & peaks_ppg_session(1:end-1)<time_end);
    %     if isempty(psqi_sg) || length(RR_seg)-length(PR_seg)>10, continue;
    %     elseif check_ppg_sqi(psqi_sg)==0, continue; end  
    %     % figure,plot(RR_seg,'*-');hold on;plot(PR_seg,'o-'); legend('HR','PR')
    %     % Analyze each segment
    %     SDNN_HR = HRV.SDNN(RR_seg,0,1); SDNN_PR = HRV.SDNN(PR_seg,0,1);
    %     RMSSD_HR = HRV.RMSSD(RR_seg,0,1); RMSSD_PR = HRV.RMSSD(PR_seg,0,1);
    %     pNN50_HR = HRV.pNN50(RR_seg,0,1); pNN50_PR = HRV.pNN50(PR_seg,0,1);
    %     ApEn_HR = HRV.ApEn(RR_seg); ApEn_PR = HRV.ApEn(PR_seg);
    % 
    %     HRV_all.Upright.RMSSD(session,sg) = RMSSD_HR; HRV_all.Upright.SDNN(session,sg) = SDNN_HR;
    %     PRV_all.Upright.RMSSD(session,sg) = RMSSD_PR; PRV_all.Upright.SDNN(session,sg) = SDNN_PR;
    %     HRV_all.Upright.pNN50(session,sg) = pNN50_HR; PRV_all.Upright.pNN50(session,sg) = pNN50_PR;
    %     HRV_all.Upright.ApEn(session,sg) = ApEn_HR; PRV_all.Upright.ApEn(session,sg) = ApEn_PR;        
    % end
    end
    % if isfield(HRV_all,'Upright')
    %     output.HRV.Upright.pNN50 = HRV_all.Upright.pNN50(:); output.PRV.Upright.pNN50 = PRV_all.Upright.pNN50(:);
    %     output.HRV.Upright.RMSSD = HRV_all.Upright.RMSSD(:); output.PRV.Upright.RMSSD = PRV_all.Upright.RMSSD(:);
    %     output.HRV.Upright.SDNN = HRV_all.Upright.SDNN(:); output.PRV.Upright.SDNN = PRV_all.Upright.SDNN(:);
    %     output.HRV.Upright.ApEn = HRV_all.Upright.ApEn(:); output.PRV.Upright.ApEn = PRV_all.Upright.ApEn(:);
    % end
    % 
    % output.HRV.Supine.pNN50 = HRV_all.Supine.pNN50(:); output.PRV.Supine.pNN50 = PRV_all.Supine.pNN50(:); 
    % output.HRV.Supine.RMSSD = HRV_all.Supine.RMSSD(:); output.PRV.Supine.RMSSD = PRV_all.Supine.RMSSD(:); 
    % output.HRV.Supine.SDNN = HRV_all.Supine.SDNN(:); output.PRV.Supine.SDNN = PRV_all.Supine.SDNN(:); 
    % output.HRV.Supine.ApEn = HRV_all.Supine.ApEn(:); output.PRV.Supine.ApEn = PRV_all.Supine.ApEn(:); 
    % 
    % save([save_path 'Results_Subject_' name],"output","id");

    % [rpc, fig, stats] = BlandAltman(nonzeros(HRV_all.Supine.SDNN(:)),nonzeros(PRV_all.Supine.SDNN(:)));
    % figure,subplot(2,2,1)
    % boxplot([HRV_all.Supine.RMSSD(:), PRV_all.Supine.RMSSD(:)]);title('RMSSD')
    % % subplot(2,2,2),plot(nonzeros(PRV_all.Supine.SDNN(:)),'*-'); hold on; plot(nonzeros(HRV_all.Supine.SDNN(:)),'*-')
    % subplot(2,2,2)
    % boxplot([HRV_all.Supine.SDNN(:),PRV_all.Supine.SDNN(:)]);title('SDNN')
    % subplot(2,2,3)
    % boxplot([HRV_all.Supine.pNN50(:),PRV_all.Supine.pNN50(:)]);title('pNN50')
    % subplot(2,2,4)
    % boxplot([HRV_all.Supine.ApEn(:),PRV_all.Supine.ApEn(:)]);title('ApEn')
    % sgtitle('Supine')
    % figure,subplot(2,2,1)
    % boxplot([HRV_all.Upright.RMSSD(:),PRV_all.Upright.RMSSD(:)]);title('RMSSD')
    % subplot(2,2,2)
    % boxplot([HRV_all.Upright.SDNN(:),PRV_all.Upright.SDNN(:)]);title('SDNN')
    % subplot(2,2,3)
    % boxplot([HRV_all.Upright.pNN50(:),PRV_all.Upright.pNN50(:)]);title('pNN50')
    % subplot(2,2,4)
    % boxplot([HRV_all.Upright.ApEn(:),PRV_all.Upright.ApEn(:)]);title('ApEn')
    % sgtitle('Upright')
    output.name = name;
    output.session_reserve = good_list;
end
function A_cleaned = delete_zeros(A)
    % delete zeros rows
    nonzero_rows = any(A, 2);% 删除全零的行
    A_cleaned = A(nonzero_rows, :);
end