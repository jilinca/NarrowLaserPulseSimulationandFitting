%% Simulation of generated gaussians with 5GHz sampling frequency
function [s_peaks, s_peaks_n, s_mp, sn_mp, one_peak_time,...                  % sampled peak values
    gen_peak, gen_ppos, geng_A, gen_FWHM,...                                % generated parameters
    n_maxval,n_ppos,...                                                     % generated parameters with noise    
    spl_mv, spl_mp, spl_el_time, all_spl_y, all_spl_x,...                   % spline fitted values
    samps, samps_n,...                                                      % All sampled values
    ypg_m,x0pg,A_pg,FWHM_pg,...                                             % parametrization fit 5 and 4GHz
    pg_mv_nr, pg_mp_nr,pg_A_nr,pg_FWHM_nr,...                               % parametrization fit noise removed 5 and 4GHz
    pg_mv2,pg_pv2,pg_A2,pg_FWHM2,...                                        % 2 GHz parmatrization
    pg_mv2_nr,pg_pv2_nr,pg_A2_nr,pg_FWHM2_nr,...                            % 2 GHz parmatrization noise removed
    pg_mv1,pg_pv1,pg_A1,pg_FWHM1,...                                        % 1 GHz parmatrization
    pg_mv1_nr,pg_pv1_nr,pg_A1_nr,pg_FWHM1_nr,...                            % 1 GHz parametrization noise removed 
    n_valm,n_vals,...                                                       % stored noise values from around the pulse (sampled data)
    pol_mv,pol_mp,pol_el_time,all_pol_y, all_pol_x,n,...                    % polynomial fit
    lm_mv,lm_mp,lm_fwhm,el_t_LM,...                                         % Levenberg-Marquardt
    yk,yk_n,xk,len_xpol,TIME_GP] =...
    fiveGHz_gen_and_samp(inter,num_p,i,sg,step)
    
n2 =nargout
%clf
%% Definition of Gaussian Parameters

% Generation of pseudo random gaussian parameters
sig1 = 0.4+(0.1)*rand(1);
% sig1 = 0.45; % meant for algolirthm developping purposes
A1 = 590+20*rand;%([590 610]);
%mu1 = randi([140 160]);

% Intervals for generation and plotting
s = 0.0125;
t =0;
%1.
xk = linspace(0.0125,inter,num_p);
%}
yk = zeros(num_p,1);                                % y-axis
% for plotting samps (can be commented away if plotting is not wanted)
sy1 = yk; sy2 = sy1; sy3 =sy1; sx =xk; % sy = zeros(num_p,3);          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_pos = find(xk>= 140 & xk <=160); % Center position
mu1 = xk(datasample(mu_pos,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generation of Gaussian
% Return pulse
yk(:,1) = A1.*exp(-(xk-mu1).^2/(2*sig1^2));
[gen_peak,p] = max(yk);
gen_ppos = xk(p);
fun_g = @(xk) A1.*exp((-(xk-mu1).^2)/(2.*sig1.^2));
% FWHM of return pulses
gen_FWHM = 2*sqrt(2*log(2))*sig1;
% Areas of return pulses
%geng_A = integral(fun_g,min(xk),max(xk));
geng_A = sqrt(2*pi)*A1*sig1;

%% Addition of gaussian noise
[yk_n, n_maxval,n_ppos] = ch_noise_modeling(yk,xk, inter,num_p);

%% SAMPLING
% Matrices for storing sampled peak amplitudes and their positions in time
s_peaks = zeros(1,sg); s_peaks_n = s_peaks;
s_mp = zeros(1,sg); sn_mp = s_mp; 

samps = zeros(num_p/step,sg);   % matrix for samples from original pulse
samps_n = samps;                % matrix for samples from noisey pulse
sp = [1 5 9 12 15];             % starting positions in yk for sampling

% SAMPLING ORIGINAL
%sampling grids (% These values are needed later for the fitting)
samps(:,1) = yk(sp(1):step:end);
samps(:,2) = yk(sp(2):step:end);
samps(:,3) = yk(sp(3):step:end);
samps(:,4) = yk(sp(4):step:end);
if sg > 4
samps(:,5) = yk(sp(5):step:end);
end
% Maximums of sampled values
[s_peaks(1,:),] = max(samps);

%SAMPLING NOISY
%sampling grids (% These values are needed later for the fitting)
samps_n(:,1) = yk_n(sp(1):step:end);
samps_n(:,2) = yk_n(sp(2):step:end);
samps_n(:,3) = yk_n(sp(3):step:end);
samps_n(:,4) = yk_n(sp(4):step:end);
if sg >4
samps_n(:,5) = yk_n(sp(5):step:end);
end

%% RETREIVAL OF MAX AMP AND POS OF SAMPLES
% Maximum sampled values and postions

[s_peaks(1,:),s_mp(1,:)]=max(samps);                        %Original
[s_peaks_n(1,:),sn_mp(1,:)] = max(samps_n);                 %Noisy
for b=1:sg
    s_mp(1,b)=(inter/num_p)*((s_mp(1,b)-1)*step+sp(b));     %Original
    sn_mp(1,b)=(inter/num_p)*((sn_mp(1,b)-1)*step+sp(b));     %Noisy
end

%Time taken to find one peak value
tic
max(samps_n(:,1));
one_peak_time =toc;

% gathering noise data from all pulses
% finding the loacation of the maximum
mpsn = find(samps_n(:,1) == max(samps_n(:,1)));
% gathering values from around the pulse for noise analysis
if step < 21 % 5 and 4GHz sampling frequency
    n_vals(:,1)...
    = vertcat(samps_n(mpsn-519:mpsn-20,1),samps_n(mpsn+20:mpsn+519,1));
    % mean of noise
    n_valm = mean(n_vals);
elseif step == 40
    n_vals(:,1)...    
    = vertcat(samps_n(mpsn-269:mpsn-10,1),samps_n(mpsn+10:mpsn+269,1));
    n_valm = mean(n_vals);
    n_vals =0;
else
    n_vals = 0;
    n_valm =20;  % This value is used in the cases of 1GHZ and 2GHZ since
    % in reality the system noise is known 
end

%% SPLINE INTERPOLATION

% Returns Maximum values, FWHM, and AREA for all sampling grids
% in_args contain both noisy waveform and orig. waveform (yk_n and yk)
% in_args contain both noisy samples and orig. samples (samps_n and samps)

[spl_mv, spl_mp, spl_el_time,all_spl_y,all_spl_x] =...
    splint2(samps,samps_n,num_p,inter,step,sg,yk,yk_n,i,xk,sp);
%tim = timeit(@()splint(samps,samps_n,num_p,inter,step,sg,yk,i))
%% POLYNOMIAL FITTING

[pol_mv,pol_mp, pol_el_time, len_xpol, all_pol_y, all_pol_x,n] =...
    poly_fit2(samps, samps_n,num_p,inter,step,sg,yk,yk_n,xk,i,sp);
%time = timeit(@()poly_fit(samps, samps_n,num_p,inter,step,sg,yk,yk_n,i,n))

%% LEVENBERG MARQUARDT FIT

[lm_mv,lm_mp,lm_fwhm,el_t_LM] = ...
    LM_fit2(samps,samps_n,yk,yk_n,xk,sg,step,inter,num_p,i,n_valm,sp);

%% Gaussian paramterization using equation (fit_gauss_two)
% this function retrieves the peak position from the samples and their
% correct locations in time, and gives them as inputs for the fit_gauss_two
% function which calculates gaussian waveform parameters
% returns mean values for multiple samples. 
[ypg_m,x0pg,A_pg,FWHM_pg,pg_mv_nr, pg_mp_nr,pg_A_nr,pg_FWHM_nr,...          % 5 and 4 GHz ret args
pg_mv2,pg_pv2,pg_A2,pg_FWHM2,pg_mv2_nr,pg_pv2_nr,pg_A2_nr,pg_FWHM2_nr,...   % 2 GHz ret args
pg_mv1,pg_pv1,pg_A1,pg_FWHM1,pg_mv1_nr,pg_pv1_nr,pg_A1_nr,pg_FWHM1_nr,...
TIME_GP]...  
= get_peaks_for_fg2(samps,samps_n,yk,yk_n,xk,sg,step,inter,num_p,i,n_valm);

end