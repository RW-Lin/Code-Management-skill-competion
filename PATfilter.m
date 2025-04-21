function [PAT, outtime] = PATfilter(PATin, time, limit)
% PATfilter: Artifact filtering for PAT sequences.
%   [PAT, time] = PATfilter(PAT, time, limit) removes artifacts using relative PAT intervals.
%   The function retains the corresponding time values for the filtered PAT.

    PAT = PATin(:); % Ensure PAT is a column vector
    outtime = time(:); % Ensure time is a column vector

    if nargin < 3
        limit = 20; % Default threshold: 20%
    end

    % Remove unreasonable values (e.g., PAT > 0.5 seconds)
    invalid_idx = PAT > 500 & PAT < 50;
    PAT(invalid_idx) = NaN;
    outtime(invalid_idx) = NaN;

    % Calculate relative differences as percentages
    pat_pct = 100 * diff(PAT) ./ PAT(1:end-1);

    % Initial artifact removal based on the specified limit
    invalid_idx = [false; abs(pat_pct) > limit];
    PAT(invalid_idx) = NaN;
    outtime(invalid_idx) = NaN;

    % Recompute relative differences after the first pass
    pat_pct = 100 * diff(PAT) ./ PAT(1:end-1);

    % Iterative filtering with decreasing thresholds
    for thr = 80:-10:limit
        invalid_idx = [false; abs(pat_pct) > thr];
        PAT(invalid_idx) = NaN;
        outtime(invalid_idx) = NaN;
        pat_pct = 100 * diff(PAT) ./ PAT(1:end-1);
    end

    % Remove NaN values from time to match PAT
    valid_idx = ~isnan(PAT);
    PAT = PAT(valid_idx);
    outtime = outtime(valid_idx);

    % figure,plot(time,PATin);hold on;plot(outtime,PAT)

end
