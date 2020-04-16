function [lm_mv,lm_mp,lm_fwhm,el_t_LM]=...
    LM_fit2(samps,samps_n,yk,yk_n,xk,sg,step,inter,num_p,c,n_valm,sp)

sf = num2str(num_p/(step*inter));
time = (inter/num_p);

% Matrices for storing parameters
lm_mv = zeros(sg,1); lm_mp = lm_mv; lm_fwhm = lm_mv; el_time_LM = lm_mv;
[slocsr,slocsc] = find(samps >= 10^-3);

pos = 1;
for i=1:sg
    nofs=sum(slocsc(:)==i);% number of sampled datapoints chosen from grid
    % devides the plotting range into intervals with the step size
    %% Parameters for plotting original generated pulse
    gpp1 = find(yk == samps(slocsr(pos,1),i));          
    if length(gpp1) >1
        gpp1 = gpp1(1);
    end
    gpp2 =find(yk == samps(slocsr(pos+nofs-1,1),i));    
    if length(gpp2) >1
        gpp2 = gpp2(end);
    end
    slocs =time*((slocsr(pos:pos+nofs-1)-1)*step+sp(i));

    xx= linspace(time*((slocsr(pos,1)-1)*step+sp(i)),...
    time*((slocsr(pos+nofs-1,1)-1)*step+sp(i)),nofs);
    xx=xx.';
    size(xx)
    yy = samps_n(slocsr(pos,1):slocsr(pos+nofs-1,1),i);
    yy = yy-n_valm;  
    % Guess function is the initial guessing function used by the LM-alg
    % Here x represents the gaussian parameters that are changed
    % The initial guesses for 'x' are given in 'x0'
    %if step < 31;
        % Iitial Guess for amplitude
        A_guess = max(yy);
        pos_guess = find(yy == A_guess);
        %Initial Guess for sigma 
        hm = A_guess/2;
        idx1 = find(yy>hm,1) +[-1 0];
        idx2 = find(yy>hm,1,'last') + [0 1];
        x1 = interp1(yy(idx1),xx(idx1),hm);
        x2 = interp1(yy(idx2),xx(idx2),hm);
        fwhm= x2-x1;
        sig_guess = fwhm/(2*sqrt(2*log(2)));
        x0 =[A_guess, xx(pos_guess,1), sig_guess];
    
    guess_fun = @(x) x(1).*exp((-(xx-x(2)).^2)./(2.*x(3).^2));
    obj=@(x) yy-guess_fun(x);
    
    lb = [550, 140, 0.3]; ub=[630,160,0.6];
    
    %options for lsqcurvefit telling which algorithm it should use
    %opt.Display ='iter';
    %opt.title = 'fitting Gaussian curve';
    opt.Jacobian = 'romberg';
    
    tic
    %call for LM-function
    x = LevenbergMarquardt(obj,x0,lb,ub,opt);
    el_time_LM(i,1) = toc;
    
    lm_mv(i,1)= x(1);
    lm_mp(i,1)= x(2);
    lm_fwhm(i,1) = 2*sqrt(2*log(2))*x(3);
    
   
    pos=pos+nofs;
end 
el_t_LM =mean(el_time_LM);
end