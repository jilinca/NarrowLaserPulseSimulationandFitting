
%% Simulations for sampling and fitting withnGHz sampling frequency
% The point of this function and sub functions is to test different fitting
% methods. 
% This main function is used to run the following function which can be used 
% either in conjuction with one another or just by them selves
% (1) fiveGHz_gen_and_samp.m: generates and samples gaussian pulse
% (2) besos.m: Adds oscillations to the end of the pulse to mimic the
% device as it functions. 
% (3) splint.m: Spline interpolation for generated pulses without
% oscillations
% (4) bslpint.m: Spline interpolation for generated pulse containg 
% oscillations 
% (5) get_peaks_for_fg2: cals fit_gauss_two.m and generates args for it
% from sampled data

function [mv_err, fwhm_err, A_err, pos_err, spl_el_time, pol_el_time,...
    len_xpol, hist_mv,hist_mp,step,sf,Time_LM,GP_time]...
    = FiveGHz_main(mv_err,fwhm_err,A_err,pos_err,n_calls)
 clf
 nargout
% Path to folder to which different figures are saved
s_path  ='PathToFolder\InWhich\FilesAreSaved\';

%% Paramters for plotting and matrices for storage of values
inter = 300;  % actual plotting interval 
n_p = 100;    % number of different pulses analyzed
%% 6000 points
%sg = 4;                     % the number of sampling grids used (max = 5)
%step = 4;                   % the sampling interval used in the grids 
%num_p = 6000; % number of points
%% 24000 points
% For all other cases except the f_s= 3GHz 24000 for 3GHz  num_p = 27000; 
num_p = 24000; 
sg = 5;                    % the number of sampling grids used (max = 5)

% the sampling interval used in the grids % 5GHz = 16, 4GHz = 20, 3GHz = 30
%, 2GHz = 40, 1GHz = 80;
step = 16;     
sf = num2str(num_p/(step*inter)); % Sampling Frequency
%(all of which are sampled with 'sg' different sampling grids)
%% MATRICES AND CELLS FOR HOLDING DATA
all_s_peaks = zeros(n_p,sg); % stores all sampled peak values.
all_sn_peaks = all_s_peaks;  % stores all sampled noise peak values.
all_s_mp = all_s_peaks;        % stores all sampled peak locations
all_sn_mp = all_s_peaks;       % stores all sampled peak locations   

all_mv_spl = all_s_peaks;    % stores all spline peak values
all_mp_spl = all_s_peaks;    % stores all spline peak value positions  
all_spl_y ={n_p,1};          % stores all spline y-values
all_spl_x ={n_p,1};          % stores all spline x-values   

all_mv_pol = all_s_peaks;    % stores all polyfit peak values
all_mp_pol = all_s_peaks;    % stores all polyfit peak value positions  
all_pol_y ={n_p,1};          % stores all polyfit y-values
all_pol_x ={n_p,1};          % stores all polyfit x-values        

gen_peaks = zeros(n_p,1);   % The maximum values of the generated pulses
gen_ppos = gen_peaks;       % Positions of max values of generated pulses
gen_A = gen_peaks;          % Areas of generated pulses
gen_FWHM = gen_peaks;       % FWHM of generated pulses

ypg_m ={n_p,1};             % peak values fit_gauss_two.m 5 and 4GHz
x0pg = ypg_m;               % peak positions, fit_gauss_two.m 5 and 4GHz
A_pg = ypg_m;               % Areas of fg2 pulses, 5 and 4GHz
FWHM_pg = ypg_m;            % FWHM of fg2 pulses, 5 and 4GHz 
pg_mv_nr= ypg_m;            % Peak fg2 values fg2 noise removed,5 and 4GHz
pg_mp_nr =ypg_m;            % Peak posoition fg2 noise removed, 5 and 4GHz
pg_A_nr=ypg_m;              % Areas of fg2 pulses noise removed, 5 and 4GHz
pg_FWHM_nr =ypg_m;          % FWHM of fg2 noise removed, 5 and 4GHz 

pg_mv2 = ypg_m;             % same as above but for 2GHz
pg_pv2 = ypg_m;
pg_FWHM2 = ypg_m;
pg_A2 =ypg_m;
pg_mv2_nr = ypg_m;
pg_pv2_nr = ypg_m;
pg_FWHM2_nr = ypg_m;
pg_A2_nr =ypg_m;

pg_mv1 = all_s_peaks;       % same as above but for 1GHz
pg_pv1 = all_s_peaks;
pg_FWHM1 = all_s_peaks;
pg_A1 = all_s_peaks;
pg_mv1_nr = all_s_peaks;
pg_pv1_nr = all_s_peaks;
pg_FWHM1_nr = all_s_peaks;
pg_A1_nr = all_s_peaks;

all_mv_lm = all_s_peaks;    % stores all LM peak values
all_mp_lm = all_s_peaks;    % stores all LM peak value positions
all_fwhm_lm = all_s_peaks;  % store all LM FWHS
el_t_LM = gen_peaks;        % Time for LM-alg to compute
    
n_vals = zeros(n_p*1000,1); % vector for noise values from all pulses (5 and 4GHz only, otherwise a zero array will be returned)
n_maxval = gen_peaks;       % maximum values of noisy pulses
n_ppos = gen_peaks;         % max position of noisy pulse

spl_el_time = gen_peaks;    % elapsed time values for spline fitting
pol_el_time = gen_peaks;    % elapsed time values for polynomial fitting
peak_el_time = gen_peaks;   % elapsed time for peak aquisition

%number of sampled points, all sampling grids for all pulses)
all_samps = zeros(num_p/step,sg*n_p);
% each pulse is stored in its own column. 
% You must give special consideration for the 3GHz case.
all_samps_n = all_samps;    % storing all noisy samples for timing purposes
len_xpol =zeros(n_p,1);     % length of polynomial x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generation of Pulses and sampling
% the function returns the max values of the generated pulses, the 
% max values of the samples, and the samples them selves
% currently returns only the samples from the last pulse being analyzed.
% This can be changed by making the return value an array similar to 
% all peaks below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GP_time = 0
temp = 0;
pos =1;
pos2 =1;
for i=1:n_p
    [all_s_peaks(i,:),all_sn_peaks(i,:),all_s_mp(i,:), all_sn_mp(i,:),...   % sampled peaks
    peak_el_time(i,1),...               
    gen_peaks(i,1),gen_ppos(i,1),gen_A(i,1),gen_FWHM(i,1),...               % generated parameters
    n_maxval(i,1),n_ppos(i,1),...                                           % generated parameters with noise
    all_mv_spl(i,:),all_mp_spl(i,:),spl_el_time(i,1),all_spl_y{i,1},...     % spline
    all_spl_x{i,1},...
    all_samps(:,pos2:pos2+sg-1),all_samps_n(:,pos2:pos2+sg-1),...           % Sampled fullwaveforms
    ypg_m{i,1}, x0pg{i,1},A_pg{i,1},FWHM_pg{i,1},...                        % parametrized
    pg_mv_nr{i,1}, pg_mp_nr{i,1},pg_A_nr{i,1},pg_FWHM_nr{i,1},...           % parametrized noise removed
    pg_mv2{i,1},pg_pv2{i,1},pg_A2{i,1},pg_FWHM2{i,1},...                    % 2 GHz parmatrization
    pg_mv2_nr{i,1},pg_pv2_nr{i,1},pg_A2_nr{i,1},pg_FWHM2_nr{i,1},...        % 2 GHz parmatrization noise removed
    pg_mv1(i,:),pg_pv1(i,:),pg_A1(i,:),pg_FWHM1(i,:),...                    % 1 GHz parmatrization
    pg_mv1_nr(i,:),pg_pv1_nr(i,:),pg_A1_nr(i,:),pg_FWHM1_nr(i,:),...        % 1 GHz parametrization noise removed 
    n_valm,n_vals(pos:pos+1000-1),...                                       % stored noise values from around the pulses (sampled data)
    all_mv_pol(i,:),all_mp_pol(i,:),pol_el_time(i,1),all_pol_y{i,1},...     % Polynomial fitting(N = DEGREE OF POLYNOMIAL)
    all_pol_x{i,1},n,...
    all_mv_lm(i,:),all_mp_lm(i,:),all_fwhm_lm(i,:),el_t_LM(i,1),...         % Levenberg Marquardt
    yk,yk_n,xk,len_xpol(i,1),GP_time]...                                            % last generated full waveform (does not store all)
        =fiveGHz_gen_and_samp(inter,num_p,i,sg,step);
    
    
    pos = pos+1000;
    pos2 = pos2+sg;
    n_valm = n_valm+temp; % storing noise values
    temp = n_valm;
    GP_time = GP_time+GP_time;
end
GP_time = GP_time/n_p;
len_xpol = mean(len_xpol);
tic
max(all_samps_n);
t_all_samps_n = toc;
%%%%%%%%%%%%%%%%%%%%%% TIMING CALC (CONTINUE THIS) %%%%%%%%%%%%%%%%%%%%%%%%
format long

mean_spl_el_time = mean(spl_el_time);
mean_pol_el_time = mean(pol_el_time);
mean_peak_el_time = mean(peak_el_time); 
comp = mean_pol_el_time/mean_spl_el_time;
% levenberg marquardt
Time_LM = mean(el_t_LM(:,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOISE REMOVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
% Evaluation of noise parameters
n_valm = n_valm/n_p;
if step < 21;
    rsm_n = sqrt(sum(abs(n_valm-n_vals).^2)/length(n_vals));
end
all_nr_peaks = all_sn_peaks-n_valm;

% Spline 
tic
all_mv_spln = all_mv_spl(:,:)-n_valm;
subs_time =toc;
  
%tic
%oneval = all_mv_spl(1,1)-n_valm;
%subs_time_2 = toc;

% Polynomial
all_mv_pol = all_mv_pol(:,:)-n_valm;

% original noisy samples are lost here 
% and noise removed samples are stored instead
for i=1:n_p
    all_spl_y{i,1} = all_spl_y{i,1}(:,:)-n_valm;
    for j=1:sg
        all_pol_y{i,1}{j,1} = all_pol_y{i,1}{j,1}(:,:)-n_valm;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PEAK VALUE ERROR CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
peak_value_errors = {5,1};

[peak_value_errors,hist_mv,p_nums] =...
    peak_err(peak_value_errors,gen_peaks,all_nr_peaks,...
    all_mv_pol,all_mv_spln,all_mv_lm, pg_mv_nr,pg_mv2_nr,...
    pg_mv1_nr,n_p,sg,step);

%Sampled max vals    
mv_err{1,1}.smp_err = mv_err{1,1}.smp_err+peak_value_errors{1,1}.smp_err;
mv_err{1,1}.smp_std = mv_err{1,1}.smp_std+peak_value_errors{1,1}.smp_std;    
mv_err{1,1}.smp_m_std = mv_err{1,1}.smp_m_std+...
    peak_value_errors{1,1}.smp_m_std;    
mv_err{1,1}.smp_p_err = mv_err{1,1}.smp_p_err+...
    peak_value_errors{1,1}.smp_p_err;
mv_err{1,1}.smp_m_RSTD = mv_err{1,1}.smp_m_RSTD+...
    peak_value_errors{1,1}.smp_m_RSTD;
%Spline max vals
mv_err{2,1}.spl_err = mv_err{2,1}.spl_err+peak_value_errors{2,1}.spl_err;
mv_err{2,1}.spl_std = mv_err{2,1}.spl_std+peak_value_errors{2,1}.spl_std;    
mv_err{2,1}.spl_m_std = mv_err{2,1}.spl_m_std+...
    peak_value_errors{2,1}.spl_m_std;    
mv_err{2,1}.spl_p_err = mv_err{2,1}.spl_p_err+...
    peak_value_errors{2,1}.spl_p_err;
mv_err{2,1}.spl_m_RSTD = mv_err{2,1}.spl_m_RSTD+...
    peak_value_errors{2,1}.spl_m_RSTD;
%Polynomial max vals
mv_err{3,1}.pol_err = mv_err{3,1}.pol_err+peak_value_errors{3,1}.pol_err;
mv_err{3,1}.pol_std = mv_err{3,1}.pol_std+peak_value_errors{3,1}.pol_std;    
mv_err{3,1}.pol_m_std = mv_err{3,1}.pol_m_std+...
    peak_value_errors{3,1}.pol_m_std;    
mv_err{3,1}.pol_p_err = mv_err{3,1}.pol_p_err+...
    peak_value_errors{3,1}.pol_p_err;
mv_err{3,1}.pol_m_RSTD = mv_err{3,1}.pol_m_RSTD+...
    peak_value_errors{3,1}.pol_m_RSTD;
%Parametrization max vals
mv_err{4,1}.pg_err = mv_err{4,1}.pg_err+peak_value_errors{4,1}.pg_err;
mv_err{4,1}.pg_std = mv_err{4,1}.pg_std+peak_value_errors{4,1}.pg_std;    
mv_err{4,1}.pg_m_std = mv_err{4,1}.pg_m_std+...
    peak_value_errors{4,1}.pg_m_std;    
mv_err{4,1}.pg_p_err = mv_err{4,1}.pg_p_err+...
    peak_value_errors{4,1}.pg_p_err;
mv_err{4,1}.pg_m_RSTD = mv_err{4,1}.pg_m_RSTD+...
    peak_value_errors{4,1}.pg_m_RSTD;
%Levnberg Marquardt errors
mv_err{5,1}.lm_err = mv_err{5,1}.lm_err+peak_value_errors{5,1}.lm_err;
mv_err{5,1}.lm_std = mv_err{5,1}.lm_std+peak_value_errors{5,1}.lm_std;    
mv_err{5,1}.lm_m_std = mv_err{5,1}.lm_m_std+...
    peak_value_errors{5,1}.lm_m_std;    
mv_err{5,1}.lm_p_err = mv_err{5,1}.lm_p_err+...
    peak_value_errors{5,1}.lm_p_err;
mv_err{5,1}.lm_m_RSTD = mv_err{5,1}.lm_m_RSTD+...
    peak_value_errors{5,1}.lm_m_RSTD;
%%}

%%%%%%%%%%%%%%%%%%%%%%%%%% FWHM CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fwhm_spl,fwhm_pol] =...
    fwhm_calc(all_spl_y,all_spl_x,all_pol_y,all_pol_x,n_p,sg);

%%%%%%%%%%%%%%%%%%%%%%% FWHM ERROR CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwhm_errors = {4,1};
[fwhm_errors] = fwhm_err_calc(fwhm_errors,gen_FWHM,fwhm_spl,fwhm_pol,...
    all_fwhm_lm,n_p,sg,p_nums);

%Spline FWHM Errors
fwhm_err{1,1}.spl_err = fwhm_err{1,1}.spl_err+fwhm_errors{1,1}.spl_err;
fwhm_err{1,1}.spl_std = fwhm_err{1,1}.spl_std+fwhm_errors{1,1}.spl_std;    
fwhm_err{1,1}.spl_m_std = fwhm_err{1,1}.spl_m_std+...
    fwhm_errors{1,1}.spl_m_std;    
fwhm_err{1,1}.spl_p_err = fwhm_err{1,1}.spl_p_err+...
    fwhm_errors{1,1}.spl_p_err;
fwhm_err{1,1}.spl_m_RSTD = fwhm_err{1,1}.spl_m_RSTD+...
    fwhm_errors{1,1}.spl_m_RSTD;
%Polynomial FWHM Errors   
fwhm_err{2,1}.pol_err = fwhm_err{2,1}.pol_err+fwhm_errors{2,1}.pol_err;
fwhm_err{2,1}.pol_std = fwhm_err{2,1}.pol_std+fwhm_errors{2,1}.pol_std;    
fwhm_err{2,1}.pol_m_std = fwhm_err{2,1}.pol_m_std+...
    fwhm_errors{2,1}.pol_m_std;    
fwhm_err{2,1}.pol_p_err = fwhm_err{2,1}.pol_p_err+...
    fwhm_errors{2,1}.pol_p_err;
fwhm_err{2,1}.pol_m_RSTD = fwhm_err{2,1}.pol_m_RSTD+...
    fwhm_errors{2,1}.pol_m_RSTD;
%Levenberg Marquardt FWHM Errors
fwhm_err{4,1}.lm_err = fwhm_err{4,1}.lm_err+fwhm_errors{4,1}.lm_err;
fwhm_err{4,1}.lm_std = fwhm_err{4,1}.lm_std+fwhm_errors{4,1}.lm_std;    
fwhm_err{4,1}.lm_m_std = fwhm_err{4,1}.lm_m_std+fwhm_errors{4,1}.lm_m_std;    
fwhm_err{4,1}.lm_p_err = fwhm_err{4,1}.lm_p_err+fwhm_errors{4,1}.lm_p_err;
fwhm_err{4,1}.lm_m_RSTD = fwhm_err{4,1}.lm_m_RSTD+...
    fwhm_errors{4,1}.lm_m_RSTD;

%% PARAMETRIZED NOISE REOMVAL AND FWHM ERR CALCULARTION
%Peak Value
for d=1:n_p
    % Parametrization
    ypg_m{d,1} = ypg_m{d,1}(:,:)-n_valm;
end

[fwhm_errors] = param_FWHM_err_calc(fwhm_errors,gen_FWHM,pg_FWHM_nr,...
    pg_FWHM2_nr,pg_FWHM1_nr,n_p,sg,step);

%Parametrization FWHM Errors
fwhm_err{3,1}.pg_err = fwhm_err{3,1}.pg_err+fwhm_errors{3,1}.pg_err;
fwhm_err{3,1}.pg_std = fwhm_err{3,1}.pg_std+fwhm_errors{3,1}.pg_std;    
fwhm_err{3,1}.pg_m_std = fwhm_err{3,1}.pg_m_std+fwhm_errors{3,1}.pg_m_std;    
fwhm_err{3,1}.pg_p_err = fwhm_err{3,1}.pg_p_err+fwhm_errors{3,1}.pg_p_err;
fwhm_err{3,1}.pg_m_RSTD = fwhm_err{3,1}.pg_m_RSTD+...
    fwhm_errors{3,1}.pg_m_RSTD;



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AREA CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Area CALCULATION
[area_spl,area_pol,area_lm] =...
    area_calc(all_mv_spl,all_mv_pol,fwhm_spl,fwhm_pol,all_spl_y,...
    all_spl_x,all_pol_y,all_pol_x,all_mv_lm, all_fwhm_lm,n_p,sg);

% Area Erreor calculation
area_errors = {4,1};
[area_errors] = ...
    area_err_calc(area_errors,gen_A,area_spl,area_pol,area_lm,pg_A_nr,pg_A2_nr,...
    pg_A1_nr,n_p,sg,step,p_nums);

%Spline 
A_err{1,1}.spl_err = A_err{1,1}.spl_err+area_errors{1,1}.spl_err;
A_err{1,1}.spl_std = A_err{1,1}.spl_std+area_errors{1,1}.spl_std;    
A_err{1,1}.spl_m_std = A_err{1,1}.spl_m_std+area_errors{1,1}.spl_m_std;    
A_err{1,1}.spl_p_err = A_err{1,1}.spl_p_err+area_errors{1,1}.spl_p_err;
A_err{1,1}.spl_m_RSTD = A_err{1,1}.spl_m_RSTD+area_errors{1,1}.spl_m_RSTD;
%Polynomial 
A_err{2,1}.pol_err = A_err{2,1}.pol_err+area_errors{2,1}.pol_err;
A_err{2,1}.pol_std = A_err{2,1}.pol_std+area_errors{2,1}.pol_std;    
A_err{2,1}.pol_m_std = A_err{2,1}.pol_m_std+area_errors{2,1}.pol_m_std;    
A_err{2,1}.pol_p_err = A_err{2,1}.pol_p_err+area_errors{2,1}.pol_p_err;
A_err{2,1}.pol_m_RSTD = A_err{2,1}.pol_m_RSTD+area_errors{2,1}.pol_m_RSTD;
%Parametrization 
A_err{3,1}.pg_err = A_err{3,1}.pg_err+area_errors{3,1}.pg_err;
A_err{3,1}.pg_std = A_err{3,1}.pg_std+area_errors{3,1}.pg_std;    
A_err{3,1}.pg_m_std = A_err{3,1}.pg_m_std+area_errors{3,1}.pg_m_std;    
A_err{3,1}.pg_p_err = A_err{3,1}.pg_p_err+area_errors{3,1}.pg_p_err;
A_err{3,1}.pg_m_RSTD = A_err{3,1}.pg_m_RSTD+area_errors{3,1}.pg_m_RSTD;
%Levenberg Marquardt
A_err{4,1}.lm_err = A_err{4,1}.lm_err+area_errors{4,1}.lm_err;
A_err{4,1}.lm_std = A_err{4,1}.lm_std+area_errors{4,1}.lm_std;    
A_err{4,1}.lm_m_std = A_err{4,1}.lm_m_std+area_errors{4,1}.lm_m_std;    
A_err{4,1}.lm_p_err = A_err{4,1}.lm_p_err+area_errors{4,1}.lm_p_err;
A_err{4,1}.lm_m_RSTD = A_err{4,1}.lm_m_RSTD+area_errors{4,1}.lm_m_RSTD;
%%}


%% Position Errors
position_errors={5,1};

[position_errors,hist_mp] = pos_err_calc(position_errors,gen_ppos,...
    all_sn_mp,all_mp_spl,all_mp_pol,all_mp_lm,pg_mp_nr,pg_pv2_nr,...
    pg_pv1_nr,n_p,sg,step);

%Sampled peak position    
pos_err{1,1}.smp_err = pos_err{1,1}.smp_err+position_errors{1,1}.smp_err;
pos_err{1,1}.smp_std = pos_err{1,1}.smp_std+position_errors{1,1}.smp_std;    
pos_err{1,1}.smp_m_std = pos_err{1,1}.smp_m_std+...
    position_errors{1,1}.smp_m_std;    
pos_err{1,1}.smp_p_err = pos_err{1,1}.smp_p_err+...
    position_errors{1,1}.smp_p_err;
pos_err{1,1}.smp_m_RSTD = pos_err{1,1}.smp_m_RSTD+...
    position_errors{1,1}.smp_m_RSTD;
%Spline peak position
pos_err{2,1}.spl_err = pos_err{2,1}.spl_err+position_errors{2,1}.spl_err;
pos_err{2,1}.spl_std = pos_err{2,1}.spl_std+position_errors{2,1}.spl_std;    
pos_err{2,1}.spl_m_std = pos_err{2,1}.spl_m_std+...
    position_errors{2,1}.spl_m_std;    
pos_err{2,1}.spl_p_err = pos_err{2,1}.spl_p_err+...
    position_errors{2,1}.spl_p_err;
pos_err{2,1}.spl_m_RSTD = pos_err{2,1}.spl_m_RSTD+...
    position_errors{2,1}.spl_m_RSTD;
%Polynomial peak position
pos_err{3,1}.pol_err = pos_err{3,1}.pol_err+position_errors{3,1}.pol_err;
pos_err{3,1}.pol_std = pos_err{3,1}.pol_std+position_errors{3,1}.pol_std;    
pos_err{3,1}.pol_m_std = pos_err{3,1}.pol_m_std+...
    position_errors{3,1}.pol_m_std;    
pos_err{3,1}.pol_p_err = pos_err{3,1}.pol_p_err+...
    position_errors{3,1}.pol_p_err;
pos_err{3,1}.pol_m_RSTD = pos_err{3,1}.pol_m_RSTD+...
    position_errors{3,1}.pol_m_RSTD;
%Parametrization peak position
pos_err{4,1}.pg_err = pos_err{4,1}.pg_err+position_errors{4,1}.pg_err;
pos_err{4,1}.pg_std = pos_err{4,1}.pg_std+position_errors{4,1}.pg_std;    
pos_err{4,1}.pg_m_std = pos_err{4,1}.pg_m_std+...
    position_errors{4,1}.pg_m_std;    
pos_err{4,1}.pg_p_err = pos_err{4,1}.pg_p_err+...
    position_errors{4,1}.pg_p_err;
pos_err{4,1}.pg_m_RSTD = pos_err{4,1}.pg_m_RSTD+...
    position_errors{4,1}.pg_m_RSTD;
%Levenberg Marquardt peak position
pos_err{5,1}.lm_err = pos_err{5,1}.lm_err+position_errors{5,1}.lm_err;
pos_err{5,1}.lm_std = pos_err{5,1}.lm_std+position_errors{5,1}.lm_std;    
pos_err{5,1}.lm_m_std = pos_err{5,1}.lm_m_std+...
    position_errors{5,1}.lm_m_std;    
pos_err{5,1}.lm_p_err = pos_err{5,1}.lm_p_err+...
    position_errors{5,1}.lm_p_err;
pos_err{5,1}.lm_m_RSTD = pos_err{5,1}.lm_m_RSTD+...
    position_errors{5,1}.lm_m_RSTD;
end
