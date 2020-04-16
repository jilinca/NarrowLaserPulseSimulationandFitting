%% THIS IS THE MAIN FUNCTION AND IT WILL CALL FiveGHz_main.m
% The purpose of this main is to retreive error results from multiple 
% runs of the first main code which performs the simulation, 
% and produce a more accurate error analysis.

s_path  ='PathToFolder\InWhich\FilesAreSaved\';
n_calls = 100;
%% MAX AMPLITUDE ERROR
mv_err={5,1}; % cell to store peak value error parameters

%struct for final sample peak value error
f_smp_err = struct('smp_err',0, 'smp_std',0, 'smp_m_std',0,...
        'smp_p_err',0,'smp_m_RSTD',0);
mv_err{1,1}=f_smp_err;
%struct for final spline errors
f_spl_err = struct('spl_err',0, 'spl_std',0, 'spl_m_std',0,...
        'spl_p_err',0,'spl_m_RSTD',0);
mv_err{2,1}=f_spl_err;
%struct for final polynomial errors
f_poly_err = struct('pol_err',0, 'pol_std',0, 'pol_m_std',0,...
        'pol_p_err',0,'pol_m_RSTD',0);
mv_err{3,1} = f_poly_err;     

%struct for final parametriztion errors
f_par_err = struct('pg_err',0, 'pg_std',0, 'pg_m_std',0,...
        'pg_p_err',0,'pg_m_RSTD',0);
mv_err{4,1} = f_par_err;

%struct for final Levenberg Marquardt errors
f_lm_err = struct('lm_err',0, 'lm_std',0, 'lm_m_std',0,...
        'lm_p_err',0,'lm_m_RSTD',0);
mv_err{5,1} = f_lm_err;

%% FWHM ERROR
%%{
fwhm_err={4,1}; % cell to store fwhm error parameters

%struct for final spline errors
ffwhm_spl_err = struct('spl_err',0, 'spl_std',0, 'spl_m_std',0,...
        'spl_p_err',0,'spl_m_RSTD',0);
fwhm_err{1,1}=ffwhm_spl_err;

%struct for final polynomial errors
ffwhm_pol_err = struct('pol_err',0, 'pol_std',0, 'pol_m_std',0,...
        'pol_p_err',0,'pol_m_RSTD',0);
fwhm_err{2,1} = ffwhm_pol_err;  

%struct for final parametrization FWHM errors
ffwhm_pg_err = struct('pg_err',0, 'pg_std',0, 'pg_m_std',0,...
        'pg_p_err',0,'pg_m_RSTD',0);
fwhm_err{3,1} = ffwhm_pg_err;
 
%struct for final Levenberg-Marquardt FWHM errors
ffwhm_lm_err = struct('lm_err',0, 'lm_std',0, 'lm_m_std',0,...
        'lm_p_err',0,'lm_m_RSTD',0);
fwhm_err{4,1} = ffwhm_lm_err;  
%%}

%% AREA ERRORS
A_err={4,1};

%struct for final spline errors
fA_spl_err = struct('spl_err',0, 'spl_std',0, 'spl_m_std',0,...
        'spl_p_err',0,'spl_m_RSTD',0);
A_err{1,1}=fA_spl_err;
%struct for final polynomial errors
fA_pol_err = struct('pol_err',0, 'pol_std',0, 'pol_m_std',0,...
        'pol_p_err',0,'pol_m_RSTD',0);
A_err{2,1} = fA_pol_err;  

%struct for final parametrization errors errors
fA_pg_err = struct('pg_err',0, 'pg_std',0, 'pg_m_std',0,...
        'pg_p_err',0,'pg_m_RSTD',0);
A_err{3,1} = fA_pg_err;

%struct for final Levenberg-Marquardt errors
fA_lm_err = struct('lm_err',0, 'lm_std',0, 'lm_m_std',0,...
        'lm_p_err',0,'lm_m_RSTD',0);
A_err{4,1} = fA_lm_err;  

%% Position Errors
pos_err={5,1};

%struct for final sample peak position error
fp_smp_err = struct('smp_err',0, 'smp_std',0, 'smp_m_std',0,...
        'smp_p_err',0,'smp_m_RSTD',0);
pos_err{1,1}=fp_smp_err;
%struct for final spline peak position errors
fp_spl_err = struct('spl_err',0, 'spl_std',0, 'spl_m_std',0,...
        'spl_p_err',0,'spl_m_RSTD',0);
pos_err{2,1}=fp_spl_err;
%struct for final polynomial peak positon errors
fp_poly_err = struct('pol_err',0, 'pol_std',0, 'pol_m_std',0,...
        'pol_p_err',0,'pol_m_RSTD',0);
pos_err{3,1} = fp_poly_err;    
%struct for final parametriztion peak position errors
fp_par_err = struct('pg_err',0, 'pg_std',0, 'pg_m_std',0,...
        'pg_p_err',0,'pg_m_RSTD',0);
pos_err{4,1} = fp_par_err;

%struct for final Levenberg_Marquardt peak positon errors
fp_lm_err = struct('lm_err',0, 'lm_std',0, 'lm_m_std',0,...
        'lm_p_err',0,'lm_m_RSTD',0);
pos_err{5,1} = fp_lm_err;

disp('Did you rember to change the spline plotting range back to original on lines 19-20 and lines: 42-43 of the splint2 function')
disp('The accuracy seems to remain the same eventhough you are plotting with 200 points instead of 24000 :)')
disp('However do more tests still')

%% Matrices for Histograms 
h_mv_smp =[];
h_mv_spl =[];
h_mv_pol =[];
h_mv_lm = [];
h_mv_pg1 =[];
h_mv_pg2 =[];
h_mv_pg3 =[];
%% Calling the simulation
for i=1:n_calls
    [mv_err,fwhm_err,A_err,pos_err,spl_el_time,pol_el_time,len_xpol,...
        hist_mv,hist_mp,step,sf,Time_LM,GP_time] =...
        FiveGHz_main(mv_err,fwhm_err,A_err,pos_err,i);
end
hist_mvf={h_mv_smp;h_mv_spl;h_mv_pol;h_mv_pg1;h_mv_pg2;h_mv_pg3};

spl_el_time =mean(spl_el_time);
pol_el_time =mean(pol_el_time);

%% Making the final divisions
% Maximum value errors
mv_err{1,1}.smp_err=mv_err{1,1}.smp_err/n_calls;
mv_err{1,1}.smp_std=mv_err{1,1}.smp_std/n_calls;
mv_err{1,1}.smp_m_std=mv_err{1,1}.smp_m_std/n_calls;
mv_err{1,1}.smp_p_err=mv_err{1,1}.smp_p_err/n_calls;
mv_err{1,1}.smp_m_RSTD=mv_err{1,1}.smp_m_RSTD/n_calls;

mv_err{2,1}.spl_err=mv_err{2,1}.spl_err/n_calls;
mv_err{2,1}.spl_std=mv_err{2,1}.spl_std/n_calls;
mv_err{2,1}.spl_m_std=mv_err{2,1}.spl_m_std/n_calls;
mv_err{2,1}.spl_p_err=mv_err{2,1}.spl_p_err/n_calls;
mv_err{2,1}.spl_m_RSTD=mv_err{2,1}.spl_m_RSTD/n_calls;

mv_err{3,1}.pol_err=mv_err{3,1}.pol_err/n_calls;
mv_err{3,1}.pol_std=mv_err{3,1}.pol_std/n_calls;
mv_err{3,1}.pol_m_std=mv_err{3,1}.pol_m_std/n_calls;
mv_err{3,1}.pol_p_err=mv_err{3,1}.pol_p_err/n_calls;
mv_err{3,1}.pol_m_RSTD=mv_err{3,1}.pol_m_RSTD/n_calls;

mv_err{4,1}.pg_err=mv_err{4,1}.pg_err/n_calls;
mv_err{4,1}.pg_std=mv_err{4,1}.pg_std/n_calls;
mv_err{4,1}.pg_m_std=mv_err{4,1}.pg_m_std/n_calls;
mv_err{4,1}.pg_p_err=mv_err{4,1}.pg_p_err/n_calls;
mv_err{4,1}.pg_m_RSTD=mv_err{4,1}.pg_m_RSTD/n_calls;

mv_err{5,1}.lm_err=mv_err{5,1}.lm_err/n_calls;
mv_err{5,1}.lm_std=mv_err{5,1}.lm_std/n_calls;
mv_err{5,1}.lm_m_std=mv_err{5,1}.lm_m_std/n_calls;
mv_err{5,1}.lm_p_err=mv_err{5,1}.lm_p_err/n_calls;
mv_err{5,1}.lm_m_RSTD=mv_err{5,1}.lm_m_RSTD/n_calls;

% FWHM Errors
fwhm_err{1,1}.spl_err=fwhm_err{1,1}.spl_err/n_calls;
fwhm_err{1,1}.spl_std=fwhm_err{1,1}.spl_std/n_calls;
fwhm_err{1,1}.spl_m_std=fwhm_err{1,1}.spl_m_std/n_calls;
fwhm_err{1,1}.spl_p_err=fwhm_err{1,1}.spl_p_err/n_calls;
fwhm_err{1,1}.spl_m_RSTD=fwhm_err{1,1}.spl_m_RSTD/n_calls;

fwhm_err{2,1}.pol_err=fwhm_err{2,1}.pol_err/n_calls;
fwhm_err{2,1}.pol_std=fwhm_err{2,1}.pol_std/n_calls;
fwhm_err{2,1}.pol_m_std=fwhm_err{2,1}.pol_m_std/n_calls;
fwhm_err{2,1}.pol_p_err=fwhm_err{2,1}.pol_p_err/n_calls;
fwhm_err{2,1}.pol_m_RSTD=fwhm_err{2,1}.pol_m_RSTD/n_calls;

fwhm_err{3,1}.pg_err=fwhm_err{3,1}.pg_err/n_calls;
fwhm_err{3,1}.pg_std=fwhm_err{3,1}.pg_std/n_calls;
fwhm_err{3,1}.pg_m_std=fwhm_err{3,1}.pg_m_std/n_calls;
fwhm_err{3,1}.pg_p_err=fwhm_err{3,1}.pg_p_err/n_calls;
fwhm_err{3,1}.pg_m_RSTD=fwhm_err{3,1}.pg_m_RSTD/n_calls;

fwhm_err{4,1}.lm_err=fwhm_err{4,1}.lm_err/n_calls;
fwhm_err{4,1}.lm_std=fwhm_err{4,1}.lm_std/n_calls;
fwhm_err{4,1}.lm_m_std=fwhm_err{4,1}.lm_m_std/n_calls;
fwhm_err{4,1}.lm_p_err=fwhm_err{4,1}.lm_p_err/n_calls;
fwhm_err{4,1}.lm_m_RSTD=fwhm_err{4,1}.lm_m_RSTD/n_calls;

% Area Errors
A_err{1,1}.spl_err=A_err{1,1}.spl_err/n_calls;
A_err{1,1}.spl_std=A_err{1,1}.spl_std/n_calls;
A_err{1,1}.spl_m_std=A_err{1,1}.spl_m_std/n_calls;
A_err{1,1}.spl_p_err=A_err{1,1}.spl_p_err/n_calls;
A_err{1,1}.spl_m_RSTD=A_err{1,1}.spl_m_RSTD/n_calls;

A_err{2,1}.pol_err=A_err{2,1}.pol_err/n_calls;
A_err{2,1}.pol_std=A_err{2,1}.pol_std/n_calls;
A_err{2,1}.pol_m_std=A_err{2,1}.pol_m_std/n_calls;
A_err{2,1}.pol_p_err=A_err{2,1}.pol_p_err/n_calls;
A_err{2,1}.pol_m_RSTD=A_err{2,1}.pol_m_RSTD/n_calls;

A_err{3,1}.pg_err=A_err{3,1}.pg_err/n_calls;
A_err{3,1}.pg_std=A_err{3,1}.pg_std/n_calls;
A_err{3,1}.pg_m_std=A_err{3,1}.pg_m_std/n_calls;
A_err{3,1}.pg_p_err=A_err{3,1}.pg_p_err/n_calls;
A_err{3,1}.pg_m_RSTD=A_err{3,1}.pg_m_RSTD/n_calls;

A_err{4,1}.lm_err=A_err{4,1}.lm_err/n_calls;
A_err{4,1}.lm_std=A_err{4,1}.lm_std/n_calls;
A_err{4,1}.lm_m_std=A_err{4,1}.lm_m_std/n_calls;
A_err{4,1}.lm_p_err=A_err{4,1}.lm_p_err/n_calls;
A_err{4,1}.lm_m_RSTD=A_err{4,1}.lm_m_RSTD/n_calls;

%Position Errors
pos_err{1,1}.smp_err=pos_err{1,1}.smp_err/n_calls;
pos_err{1,1}.smp_std=pos_err{1,1}.smp_std/n_calls;
pos_err{1,1}.smp_m_std=pos_err{1,1}.smp_m_std/n_calls;
pos_err{1,1}.smp_p_err=pos_err{1,1}.smp_p_err/n_calls;
pos_err{1,1}.smp_m_RSTD=pos_err{1,1}.smp_m_RSTD/n_calls;

pos_err{2,1}.spl_err=pos_err{2,1}.spl_err/n_calls;
pos_err{2,1}.spl_std=pos_err{2,1}.spl_std/n_calls;
pos_err{2,1}.spl_m_std=pos_err{2,1}.spl_m_std/n_calls;
pos_err{2,1}.spl_p_err=pos_err{2,1}.spl_p_err/n_calls;
pos_err{2,1}.spl_m_RSTD=pos_err{2,1}.spl_m_RSTD/n_calls;

pos_err{3,1}.pol_err=pos_err{3,1}.pol_err/n_calls;
pos_err{3,1}.pol_std=pos_err{3,1}.pol_std/n_calls;
pos_err{3,1}.pol_m_std=pos_err{3,1}.pol_m_std/n_calls;
pos_err{3,1}.pol_p_err=pos_err{3,1}.pol_p_err/n_calls;
pos_err{3,1}.pol_m_RSTD=pos_err{3,1}.pol_m_RSTD/n_calls;

pos_err{4,1}.pg_err=pos_err{4,1}.pg_err/n_calls;
pos_err{4,1}.pg_std=pos_err{4,1}.pg_std/n_calls;
pos_err{4,1}.pg_m_std=pos_err{4,1}.pg_m_std/n_calls;
pos_err{4,1}.pg_p_err=pos_err{4,1}.pg_p_err/n_calls;
pos_err{4,1}.pg_m_RSTD=pos_err{4,1}.pg_m_RSTD/n_calls;

pos_err{5,1}.lm_err=pos_err{5,1}.lm_err/n_calls;
pos_err{5,1}.lm_std=pos_err{5,1}.lm_std/n_calls;
pos_err{5,1}.lm_m_std=pos_err{5,1}.lm_m_std/n_calls;
pos_err{5,1}.lm_p_err=pos_err{5,1}.lm_p_err/n_calls;
pos_err{5,1}.lm_m_RSTD=pos_err{5,1}.lm_m_RSTD/n_calls;


