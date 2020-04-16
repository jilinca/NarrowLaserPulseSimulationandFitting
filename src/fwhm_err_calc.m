function [fwhm_errors]=...
    fwhm_err_calc(fwhm_errors,gen_FWHM,fwhm_spl,fwhm_pol,fwhm_lm,n_p,sg,p_nums)
    %% SPLINE
    % Errors and std for Retreived FWHM 
    spl_errors = struct('spl_err',0, 'spl_std',0, 'spl_m_std',0,...
    'spl_p_err',0,'spl_m_RSTD',0);
    %mean error
    dif_spl={n_p,1};                       
    m_dif_spl =zeros(n_p,1);
    %percentage mean error
    p_dif_spl = {n_p,1};                   
    pm_dif_spl =zeros(n_p,1);
    %standard deviation
    s_spl=zeros(n_p,1);
    for i=1:n_p
       %mean error
       dif_spl{i,1} = abs(gen_FWHM(i,1)-fwhm_spl(i,:));    
       m_dif_spl(i,1) = mean(dif_spl{i,1}(1,:));                             % mean differences for each pulse from peak value  
       %percentage mean error
       p_dif_spl{i,1} = (abs(gen_FWHM(i,1)-fwhm_spl(i,:))...
           +gen_FWHM(i,1))/gen_FWHM(i,1) -1;
       pm_dif_spl(i,1) = mean(p_dif_spl{i,1}(1,:));
       % Standard deviation from true value of zero differnece
       s_spl(i,1) = sum((dif_spl{i,1}(1,:)).^2);
    end
    %Mean Error
    mf_spl =  mean(m_dif_spl(:,1));  
    spl_errors.spl_err = mf_spl;                   % mean difference from zero error 
    % Percentual mean error
    spl_errors.spl_p_err = mean(pm_dif_spl(:,1)); 
    % Standard deviation from true value of zero difference
    spl_errors.spl_std = sqrt(sum(s_spl)/(sg*n_p -1));
    % Standard deviation from mean differences
    for j=1:n_p
        %s_spl(i,1) = sum((m_dif_spl(i,1)-dif_spl{i,1}(1,:)).^2);
        s_spl(i,1) = sum((mf_spl-dif_spl{i,1}(1,:)).^2);
    end
    std_spl =sqrt(sum(s_spl)/(sg*n_p-1));
    %RTSD
    RSTD_spl = (std_spl/mf_spl)*100;
    spl_errors.spl_m_std = sqrt(sum(s_spl)/(sg*n_p -1));
    % Percentage or relative standard deviation (RSTD)
    spl_errors.spl_m_RSTD = RSTD_spl;
    fwhm_errors{1,1} = spl_errors;
    
    %% POLYNOMIAL
    % Errors and std for Retreived FWHM 
    pol_errors = struct('pol_err',0, 'pol_std',0, 'pol_m_std',0,...
    'pol_p_err',0,'pol_m_RSTD',0);
    
    %mean error
    dif_pol={n_p,1};                       
    m_dif_pol =zeros(n_p,1);
    %percentage mean error
    p_dif_pol = {n_p,1};                   
    pm_dif_pol =zeros(n_p,1);
    %standard deviation
    s_pol=zeros(n_p,1);
    
    for i=1:n_p
       %mean error
       dif_pol{i,1} = abs(gen_FWHM(i,1)-fwhm_pol(i,:));    
       m_dif_pol(i,1) = mean(dif_pol{i,1}(1,:));                             % mean differences for each pulse from peak value  
       %percentage mean error
       p_dif_pol{i,1} = (abs(gen_FWHM(i,1)-fwhm_pol(i,:))...
           +gen_FWHM(i,1))/gen_FWHM(i,1) -1;
       pm_dif_pol(i,1) = mean(p_dif_pol{i,1}(1,:));
       % Standard deviation from true value of zero differnece
       s_pol(i,1) = sum((dif_pol{i,1}(1,:)).^2);
    end

    %Mean Error
    mf_pol =  mean(m_dif_pol(:,1));  
    pol_errors.pol_err = mf_pol;                   % mean difference from zero error 
    % Percentual mean error
    pol_errors.pol_p_err = mean(pm_dif_pol(:,1)); 
    % Standard deviation from true value of zero difference
    pol_errors.pol_std = sqrt(sum(s_pol)/(sg*n_p -1));
    % Standard deviation from mean differences
    for j=1:n_p
        %s_pol(i,1) = sum((m_dif_pol(i,1)-dif_pol{i,1}(1,:)).^2);
        s_pol(i,1) = sum((mf_pol-dif_pol{i,1}(1,:)).^2);
    end
    std_pol =sqrt(sum(s_pol)/(sg*n_p-1));
    %RTSD
    RSTD_pol = (std_pol/mf_pol)*100;
    
    pol_errors.pol_m_std = sqrt(sum(s_pol)/(sg*n_p -1));
    % Percentage or relative standard deviation (RSTD)
    pol_errors.pol_m_RSTD = RSTD_pol;
    
    fwhm_errors{2,1} = pol_errors;
    
    %% LEVENBERG MARQUARDT
    % Errors and std for Retreived FWHM 
    lm_errors = struct('lm_err',0, 'lm_std',0, 'lm_m_std',0,...
    'lm_p_err',0,'lm_m_RSTD',0);
    
    l_nums=length(p_nums);
    %mean error
    dif_lm={n_p,1};                       
    m_dif_lm =zeros(n_p,1);
    %percentage mean error
    p_dif_lm = {n_p,1};                   
    pm_dif_lm =zeros(n_p,1);
    %standard deviation
    s_lm=zeros(n_p,1);
    % Removing erroneus fitting data from LM-algorithm
    gen2_FWHM = [];
    fwhm2_lm = [];
    for k =1:l_nums
        gen2_FWHM = [gen2_FWHM;gen_FWHM(p_nums(k),1)];
        fwhm2_lm = [fwhm2_lm; fwhm_lm(p_nums(k),1)];
    end
    
    for i=1:l_nums
       %mean error
       dif_lm{i,1} = abs(gen2_FWHM(i,1)-fwhm2_lm(i,:));    
       m_dif_lm(i,1) = mean(dif_lm{i,1}(1,:));                             % mean differences for each pulse from peak value  
       %percentage mean error
       p_dif_lm{i,1} = (abs(gen2_FWHM(i,1)-fwhm2_lm(i,:))...
           +gen2_FWHM(i,1))/gen2_FWHM(i,1) -1;
       pm_dif_lm(i,1) = mean(p_dif_lm{i,1}(1,:));
       % Standard deviation from true value of zero differnece
       s_lm(i,1) = sum((dif_lm{i,1}(1,:)).^2);
    end

    %Mean Error
    mf_lm =  mean(m_dif_lm(:,1));  
    lm_errors.lm_err = mf_lm;                   % mean difference from zero error 
    % Percentual mean error
    lm_errors.lm_p_err = mean(pm_dif_lm(:,1)); 
    % Standard deviation from true value of zero difference
    lm_errors.lm_std = sqrt(sum(s_lm)/(sg*l_nums -1));
    % Standard deviation from mean differences
    for j=1:l_nums
        %s_lm(i,1) = sum((m_dif_lm(i,1)-dif_lm{i,1}(1,:)).^2);
        s_lm(i,1) = sum((mf_lm-dif_lm{i,1}(1,:)).^2);
    end
    std_lm =sqrt(sum(s_lm)/(sg*l_nums-1));
    %RTSD
    RSTD_lm = (std_lm/mf_lm)*100;
    
    lm_errors.lm_m_std = sqrt(sum(s_lm)/(sg*l_nums -1));
    % Percentage or relative standard deviation (RSTD)
    lm_errors.lm_m_RSTD = RSTD_lm;
    
    fwhm_errors{4,1} = lm_errors;
    

end