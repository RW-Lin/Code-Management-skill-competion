function [output,good_ratio] = check_ppg_sqi(psqi)
    psqi(isnan(psqi)) = 0;
    if mean(psqi)>.8, output = 1;
        idx= kmeans(psqi(:),2); ratio = sum(idx==2)/length(idx);
        if ratio > 0.5, good_ratio = ratio;
        else, good_ratio = 1-ratio;end
    % elseif mean(psqi)<0.5,output = 0;good_ratio=0;
    else
        idx= kmeans(psqi(:),2); ratio = sum(idx==2)/length(idx);
        if ratio > 0.8, output = 1; good_ratio = ratio;
        elseif ratio < 0.2, output = 1; good_ratio = 1-ratio;
        else, output = 0; good_ratio=0;
        end
    end
end