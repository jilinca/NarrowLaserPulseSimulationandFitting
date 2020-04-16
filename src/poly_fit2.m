function [yy_m,xx_m,pol_el_time,len_xpol,all_y,all_x,n] = ...
    poly_fit2(samps, samps_n,num_p,inter,step,sg,yk,yk_n,xk,c,sp)

len_xpol =0;
fit_points = 1000;
sf = num2str(num_p/(step*inter));
yy_m = zeros(sg,1); xx_m = yy_m; elTime = yy_m;
all_y ={sg,1}; all_x ={sg,1};
time = inter/num_p; % defines the plotting range correctly
[slocsr,slocsc] = find(samps >= 10);
pos = 1;
for i=1:sg
    nofs =sum(slocsc(:)==i);% number of sampled datapoints chosen from grid
    n = nofs-1;
    %% Parameters for plotting original generated pulse
    % generated pulse starting position
    gpp1 = find(yk == samps(slocsr(pos,1),i)); 
    if length(gpp1) >1
        gpp1 = gpp1(1);
    end
    % generated pulse ending position
    gpp2 =find(yk == samps(slocsr(pos+nofs-1,1),i)); 
    if length(gpp2) >1
        gpp2 = gpp2(end);
    end    
    %% Just testing
    slocs =time*((slocsr(pos:pos+nofs-1)-1)*step+sp(i));
    xax = linspace(slocs(1,1),slocs(length(slocs),1),fit_points);
    all_x{i,1} = xax;
    len=length(xax);
    len_xpol= len_xpol+len;
    xax3 = linspace(slocs(1,1),slocs(length(slocs),1),length(slocs));    
    
    tic
    [ply,~,mu] = polyfit(slocs, samps_n(slocsr(pos,1):...
        slocsr(pos+nofs-1,1),i),n);
    plyy = polyval(ply,xax,[],mu);
    elTime(i,1) = toc;
    [yy_m(i,1),p] = max(plyy);
    xx_m(i,1) = xax(p);
    all_y{i,1}= plyy;

    %Postion increment to analyze correct set of samples
    pos = pos+nofs; 
end
pol_el_time = mean(elTime);
len_xpol = len_xpol/sg;
end