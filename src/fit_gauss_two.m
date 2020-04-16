%% Fitting a gaussian using an analytical solution
function [y0,x0,sigma,peak_volume,FWHM] = fit_gauss_two(peak_pos, peak_vals)
    % Array for Full waveform 

    i=2;
    % extract three (x,y) points from the top of the peak
    x1 = peak_pos(i-1);
    lny1 = log(peak_vals(i-1));
    x2 = peak_pos(i);
    y2 = peak_vals(i);  in comaprison to ParametrizeWaveform2015
    lny2 = log(y2);
    x3 = peak_pos(i+1);
    lny3 = log(peak_vals(i+1));
    % Parametrization equations
    x0 = 0.5*((x1^2*(lny3-lny2)+x2^2*(lny1-lny3)+x3^2*(lny2-lny1))...
         /(x1*(lny3-lny2)+x2*(lny1-lny3)+x3*(lny2-lny1)));
    sigma = sqrt(0.5 * ((x2-x0)^2 - (x1-x0)^2) / (lny1-lny2) );
    y0 = y2 * exp((x2-x0)^2 / (2*sigma^2));
    FWHM = 2*sqrt(2*log(2))*sigma;
    peak_volume = y0*sqrt(2*pi)*sigma;  
end