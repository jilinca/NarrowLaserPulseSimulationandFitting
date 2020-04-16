%% Add noise to waveform
% Generates gaussian distributed white noise
function [waveform, n_maxval, n_ppos] = ch_noise_modeling(waveform,xk,...
    inter,num_p)
% gaussian whight noise
sigma = 2;
mu = 20;
noise = sigma*randn(length(waveform),1)+mu;
waveform = waveform(:,1)+noise(:,1);
[n_maxval,n_ppos] = max(waveform);
% The first maximum is chosen if more than one exists 
if length(n_ppos)>1
    n_ppos = n_ppos(1);
    n_maxval = waveform(n_ppos);
end    
%Final maximum position scaled to correct time
n_ppos =xk(n_ppos);
end