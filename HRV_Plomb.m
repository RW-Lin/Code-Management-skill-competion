function [LFpower,HFpower,LFHFratio] = HRV_Plomb(RR_interval)


[freq] = plomb(RR_interval);
LFpower = abs(freq<.15 && freq>0.04).^2; HFpower = abs(freq>.15 && freq<0.4).^2;

LFHFratio = LFpower./HF_power;
end