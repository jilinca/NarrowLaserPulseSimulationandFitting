function [yy_m, xx_m,el_m_time,all_yy_n,all_xx] =...
    splint2(samps,samps_n,num_p,inter,step,sg,yk,yk_n,c,xk,sp)
% This function fits a spline to the sampled data and scales the 
% samples to the correct locations in time. The fitted wavefoms are passed
% to FiveGHz_main and the noise is removed. After this the FWHM and Area
% are calculated. Using the calc_FWHM function.
fit_points =1000;
%Save Path for images
s_path  = '\\tsclient\C\Users\jil\Desktop\Kaivos\Thesis\Essay_and_Img\';
sf = num2str(num_p/(step*inter));
time = (inter/num_p);
yy_m = zeros(sg,1); xx_m = yy_m; elTime = yy_m;
% Matrices for storing all fitted pulse values 
all_yy_n = zeros(fit_points,sg); all_xx = all_yy_n;
% Find points in the sampled values that are higher than a threshold
% for orig gen pulse.
[slocsr,slocsc] = find(samps >= 10^-3);
pos = 1;
for i=1:sg
    nofs=sum(slocsc(:)==i);% number of sampled datapoints chosen from grid 
   
    %% Parameters for plotting original generated pulse
    gpp1 = find(yk == samps(slocsr(pos,1),i));           
    if length(gpp1) >1
        gpp1 = gpp1(1);
    end
    gpp2 =find(yk == samps(slocsr(pos+nofs-1,1),i));    
    if length(gpp2) >1
        gpp2 = gpp2(end);
    end   
    %% The spline fittig procedure
    % This is for plotting the spline
    xx= linspace(time*((slocsr(pos,1)-1)*step+sp(i)),time*((slocsr(pos+nofs-1,1)-1)*step+sp(i)),fit_points);
    slocs =time*((slocsr(pos:pos+nofs-1)-1)*step+sp(i));
    % spline fitting for orig. gen pulse
    %yy = spline(slocs(:),samps(slocsr(pos,1):slocsr(pos+nofs-1,1),i),xx);  
    % spline fitting for noisy orig. pulse
    tic
    yy_n=spline(slocs(:),samps_n(slocsr(pos,1):slocsr(pos+nofs-1,1),i),xx);
    elTime(i,1) = toc;
    all_yy_n(:,i) = yy_n;
    all_xx(:,i) = xx;
    pos = pos+nofs;            % Indexing the array at the correct position
    %[yy_m(i,1), xx_m(i,1)] = max(yy);     % Stroing maximum values (orig.)
    [yy_m(i,1),p] = max(yy_n);             % Storing max values (noise)
    xx_m(i,1) = xx(p);
end
el_m_time = mean(elTime);
end

