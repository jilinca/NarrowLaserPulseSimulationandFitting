%% Get peak values to call fit_gauss 2
% gets three different peak values for all sampling grids for each pulse.
% the peak values for eachset of three will have the same peak value chosen
% and the two surrounding values will change and will be at equal distances
% on both sides of the sampled peak value.
% each different sampling grid will produce a different peak value as the
% heighest value in the peak_vals matrix
% c keeps track of which pulse is being processeed
function [ypg_m, x0pg, A_pg, FWHM_pg,pg_mv_nr, pg_mp_nr,pg_A_nr,...
    pg_FWHM_nr,pg_mv2,pg_pv2,pg_A2,pg_FWHM2,pg_mv2_nr,pg_pv2_nr,...
    pg_A2_nr,pg_FWHM2_nr,pg_mv1,pg_pv1,pg_A1,pg_FWHM1,pg_mv1_nr,...
    pg_pv1_nr,pg_A1_nr,pg_FWHM1_nr, TIME_GP]...    
    =get_peaks_for_fg2(samps,samps_n,yk,yk_n,xk,sg,step,inter, num_p,...
    c,n_valm)
sf = num2str(num_p/(step*inter));
s_path  = '\\tsclient\C\Users\jil\Desktop\Kaivos\Thesis\Essay_and_Img\';
%% Matrix initalization
% each row of the matrices below will hold the respective gaussian
% parameters (y0, x0 and standard deviation (sig)) calculated from
% different sampling points of the same pulse. there will be five rows
% all together due to sg=5 different sampling grids.
y0_fg2 =zeros(sg ,3); x0_fg2 = y0_fg2; sig_fg2 = y0_fg2;
% matrix for storing maximum parameteried y values( pg stands for
GP_time =0;
cnt =0;
% Retreived values using noisy data
ypg_m = zeros(sg,3); x0pg= ypg_m; FWHM_pg=ypg_m; A_pg = ypg_m;              % 5 sampling grids 3 different samples chosen (5 and 4 GHz)
pg_mv2 = zeros(sg,2); pg_pv2 = pg_mv2; pg_FWHM2=pg_mv2; pg_A2 = pg_mv2;     % 5 sampling grids 2 point grids (2GHz)
pg_mv1 = zeros(1,sg); pg_pv1 = pg_mv1; pg_FWHM1=pg_mv1; pg_A1 = pg_mv1;     % 5 sampling grids 1 point grid (1GHz)
% Retreived values using noise removed data
pg_mv_nr =ypg_m; pg_mp_nr =ypg_m; pg_A_nr = ypg_m; pg_FWHM_nr = ypg_m;      % 5 sampling grids 3 different samples chosen (5 and 4 GHz)
pg_mv2_nr = zeros(sg,2); pg_pv2_nr = pg_mv2;                                % 5 sampling grids 2 point grids (2GHz)
pg_FWHM2_nr=pg_mv2; pg_A2_nr = pg_mv2;                                      % 5 sampling grids 2 point grids (2GHz)
pg_mv1_nr = pg_mv1; pg_pv1_nr = pg_mv1;                                     % 5 sampling grids 1 point grid (1GHz)
pg_FWHM1_nr=pg_mv1; pg_A1_nr = pg_mv1;                                      % 5 sampling grids 1 point grid (1GHz)
%% Finding the Peak values
[slocsr,slocsc] = find(samps >= 10^-3);
pos = 1;
peak_vals = zeros(3,3); peak_pos = peak_vals;
% arrays for the 2GHz sampling frequncey case
pv_2GHZ = zeros(3,2); pp_2GHZ = zeros(2,3);
% arrays for the 1GHz sampling frequncey case
pv_1GHZ = zeros(3,1); pp_1GHZ = zeros(1,3);

for i=1:sg
    nofs =  sum(slocsc(:)==i);
    % Cuts part of sampled values corresponding to slocs_r 
    %cut_s = samps(slocsr(pos):slocsr(pos+nofs-1),i); % ORIGINAL 
    cut_s = samps_n(slocsr(pos):slocsr(pos+nofs-1),i); % NOISY
    
    % find peak values from different distances
    [mv,mps] = max(cut_s);          % maximum value and its position
    peak_vals(2,:) = mv;            % sampling frequencies of 4 and 5GHz    
    pv_2GHZ(2,:) = mv;               % 2GHz sampling frequency 
    pv_1GHZ(2,1) = mv;              % 1GH sampling frequency
    s_dists1 = [1 3 5];             % sampling frequencies of 4 and 5GHz
    s_dists2 = [1 2];               % sampling frequncies of 2GHz  
    %% Forming the point grids for different sampling frequencies 
    if step < 21
        %closest to peak
        peak_vals(1,1) = cut_s(mps-s_dists1(1));
        peak_vals(3,1) = cut_s(mps+s_dists1(1));
        % middle
        peak_vals(1,2) = cut_s(mps-s_dists1(2));
        peak_vals(3,2) = cut_s(mps+s_dists1(2));
        % farthest from peak
        peak_vals(1,3) = cut_s(mps-s_dists1(3));
        peak_vals(3,3) = cut_s(mps+s_dists1(3));
    elseif step == 40 || step == 30;
    %elseif step == 40;
        % closest to peak
        pv_2GHZ(1,1) = cut_s(mps-s_dists2(1));
        pv_2GHZ(3,1) = cut_s(mps+s_dists2(1));
        % farthest from peak
        pv_2GHZ(1,2) = cut_s(mps-s_dists2(2));
        pv_2GHZ(3,2) = cut_s(mps+s_dists2(2));
    elseif step > 75;
        pv_1GHZ(1,1) = cut_s(mps-1);
        pv_1GHZ(3,1) = cut_s(mps+1);
    end
    peak_vals_nr = peak_vals -20; % Average noise removed from peak values
    pv_2GHZ_nr = pv_2GHZ-20;
    pv_1GHZ_nr = pv_1GHZ-20;
    
    %% Parameters for finding correct position for peak_vals
    gpp1 = find(yk == samps(slocsr(pos,1),i)); 
    if length(gpp1) >1
        gpp1 = gpp1(1);
    end
    gpp2 =find(yk == samps(slocsr(pos+nofs-1,1),i)); 
    if length(gpp2)>1
        gpp2 = gpp2(end);
    end
    % ORIGINAL
    %loc = find(yk == peak_vals(2,1)); 
    % NOISY
    loc = find(yk_n == peak_vals(2,1)); 
    if length(loc) >1
        loc = loc(1);
    end  
    % finds the positions of the chosen sampled points in time and
    % discrimintes between sampling frequency 
    peak_pos(:,2) = xk(loc);    % 5 and 4GHz sampling frequency
    pp_2GHZ(:,2) = xk(loc);     % 2 GHz sampling frequency  
    if step < 21;
        for r=1:length(peak_pos(1,:))
            peak_pos(r,1) = xk(loc-s_dists1(r)*step);
            peak_pos(r,3) = xk(loc+s_dists1(r)*step);
        end
    elseif step == 40 || step == 30
    %elseif step == 40;        
        for r=1:length(pp_2GHZ(:,1))
            pp_2GHZ(r,1) = xk(loc-s_dists2(r)*step);
            pp_2GHZ(r,3) = xk(loc+s_dists2(r)*step);    
        end
    elseif step >75
        pp_1GHZ(1,1) = xk(loc-step);
        pp_1GHZ(1,2) = xk(loc); 
        pp_1GHZ(1,3) = xk(loc+step);
    end 
    % find peak postion values in time
    start = xk(gpp1); 
    finish = xk(gpp2);
    % defines calculation range for paraetrized pulse. 
    xs = start:(finish-start)/(gpp2-gpp1):finish;
    %% Calls for the parametrization function
    if step < 21
        for  r=1:length(peak_pos(1,:))
            [y0,x0,sig,peak_volume,FWHM]...
                = fit_gauss_two(peak_pos(r,:), peak_vals(:,r));
            ypg = y0.*exp((-(xs-x0).^2)/(2.*sig.^2));
            [ypg_m(i,r)] = max(ypg); % return value for peak
            x0pg(i,r) = x0;
            FWHM_pg(i,r) = FWHM;
            A_pg(i,r) = peak_volume;
            % Recalculation with noise removed from peak vals
            tic
            [y0,x0,sig,pg_A_nr(i,r),pg_FWHM_nr(i,r)]...
                = fit_gauss_two(peak_pos(r,:), peak_vals_nr(:,r));
            ypg_nr = y0.*exp((-(xs-x0).^2)/(2.*sig.^2));
            GP_el_time =toc;
            GP_time = GP_time+GP_el_time;
            cnt =cnt+1;
            
            pg_mv_nr(i,r) =max(ypg_nr);
            pg_mp_nr(i,r) =x0;
        end
    elseif step == 40 || step == 30
         for  r=1:length(pp_2GHZ(:,1))
        % REMEBER TO TIME "FIT_GAUSS_TWO
            [y0,x0,sig,pg_A2(i,r),pg_FWHM2(i,r)]...
                = fit_gauss_two(pp_2GHZ(r,:), pv_2GHZ(:,r));
            ypg = y0.*exp((-(xs-x0).^2)/(2.*sig.^2)); 
            [pg_mv2(i,r)] = max(ypg); % return value for peak
            pg_pv2(i,r) = x0;
            % Recalculation with noise removed from peak vals
            tic
            [y0,x0,sig,pg_A2_nr(i,r),pg_FWHM2_nr(i,r)]...
                = fit_gauss_two(pp_2GHZ(r,:), pv_2GHZ_nr(:,r));
            ypg_nr = y0.*exp((-(xs-x0).^2)/(2.*sig.^2));
            pg_mv2_nr(i,r) =max(ypg_nr);
            pg_pv2_nr(i,r) =x0;
            GP_el_time =toc;
            GP_time = GP_time+GP_el_time;
            cnt =cnt+1;
         end
    elseif step > 75
        [y0,x0,sig,pg_A1(1,i),pg_FWHM1(1,i)]= fit_gauss_two(pp_1GHZ, pv_1GHZ);
        ypg = y0.*exp((-(xs-x0).^2)/(2.*sig.^2)); 
        pg_mv1(1,i) = max(ypg); % return value for peak
        pg_pv1(1,i) = x0;
        % Recalculation with noise removed from peak vals
        tic
        [y0,x0,sig,pg_A1_nr(1,i),pg_FWHM1_nr(1,i)]...
                = fit_gauss_two(pp_1GHZ, pv_1GHZ_nr);
        ypg_nr = y0.*exp((-(xs-x0).^2)/(2.*sig.^2));
                    GP_el_time =toc;
        GP_time = GP_time+GP_el_time;
        cnt =cnt+1;
        pg_mv1_nr(1,i) =max(ypg_nr);
        pg_pv1_nr(1,i) =x0;
    end
    pos = pos+nofs;
end
TIME_GP = GP_time/cnt;
end