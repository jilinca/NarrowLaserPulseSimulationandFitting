function [area_errors] =...
    area_err_calc(area_errors,gen_A,area_spl,area_pol,area_lm,...
    pg_A_nr,pg_A2_nr,pg_A1_nr, n_p,sg,step,p_nums)
    
    %% Errors and std for Retreived Areas
    %SPLINE
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
       dif_spl{i,1} = abs(gen_A(i,1)-area_spl(i,:));    
       m_dif_spl(i,1) = mean(dif_spl{i,1}(1,:));                             % mean differences for each pulse from peak value  
       %percentage mean error
       p_dif_spl{i,1} = (abs(gen_A(i,1)-area_spl(i,:))...
           +gen_A(i,1))/gen_A(i,1) -1;
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
        s_spl(i,1) = sum((mf_spl-dif_spl{i,1}(1,:)).^2);
    end
    std_spl =sqrt(sum(s_spl)/(sg*n_p-1));
    %RTSD
    RSTD_spl = (std_spl/mf_spl)*100;
    
    spl_errors.spl_m_std = sqrt(sum(s_spl)/(sg*n_p -1));
    % Percentage or relative standard deviation (RSTD)
    spl_errors.spl_m_RSTD = RSTD_spl;
    area_errors{1,1} = spl_errors;
    %% POLYNOMIAL
    % Errors and std for Retreived AREA 
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
       dif_pol{i,1} = abs(gen_A(i,1)-area_pol(i,:));    
       m_dif_pol(i,1) = mean(dif_pol{i,1}(1,:));                             % mean differences for each pulse from peak value  
       %percentage mean error
       p_dif_pol{i,1} = (abs(gen_A(i,1)-area_pol(i,:))...
           +gen_A(i,1))/gen_A(i,1) -1;
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
        s_pol(i,1) = sum((m_dif_pol(i,1)-dif_pol{i,1}(1,:)).^2);
        s_pol(i,1) = sum((mf_pol-dif_pol{i,1}(1,:)).^2);
    end
    std_pol =sqrt(sum(s_pol)/(sg*n_p-1));
    %RTSD
    RSTD_pol = (std_pol/mf_pol)*100;
    pol_errors.pol_m_std = sqrt(sum(s_pol)/(sg*n_p -1));
    % Percentage or relative standard deviation (RSTD)
    pol_errors.pol_m_RSTD = RSTD_pol;
    area_errors{2,1} = pol_errors;
    %% LEVENBERG MARQUARDT
    % Errors and std for Retreived AREA 
    lm_errors = struct('lm_err',0, 'lm_std',0, 'lm_m_std',0,...
    'lm_p_err',0,'lm_m_RSTD',0);
    %mean error
    dif_lm={n_p,1};                       
    m_dif_lm =zeros(n_p,1);
    %percentage mean error
    p_dif_lm = {n_p,1};                   
    pm_dif_lm =zeros(n_p,1);
    %standard deviation
    s_lm=zeros(n_p,1);
    %Removing pulse data where Levenberg Marquardt did not succeed.
    l_nums = length(p_nums(1,:));
    gen2_A = [];
    area2_lm = [];
    for k =1:l_nums
        gen2_A = [gen2_A;gen_A(p_nums(k),1)];
        area2_lm = [area2_lm; area_lm(p_nums(k),1)];
    end
    for i=1:l_nums
       %mean error
       dif_lm{i,1} = abs(gen2_A(i,1)-area2_lm(i,:));    
       m_dif_lm(i,1) = mean(dif_lm{i,1}(1,:));                             % mean differences for each pulse from peak value  
       %percentage mean error
       p_dif_lm{i,1} = (abs(gen2_A(i,1)-area2_lm(i,:))...
           +gen2_A(i,1))/gen2_A(i,1) -1;
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
    lm_errors.lm_std = sqrt(sum(s_lm)/(sg*n_p -1));
    % Standard deviation from mean differences
    for j=1:l_nums
        %s_lm(i,1) = sum((m_dif_lm(i,1)-dif_lm{i,1}(1,:)).^2);
        s_lm(i,1) = sum((mf_lm-dif_lm{i,1}(1,:)).^2);
    end
    std_lm =sqrt(sum(s_lm)/(sg*n_p-1));
    %RTSD
    RSTD_lm = (std_lm/mf_lm)*100;
    lm_errors.lm_m_std = sqrt(sum(s_lm)/(sg*n_p -1));
    % Percentage or relative standard deviation (RSTD)
    lm_errors.lm_m_RSTD = RSTD_lm;
    area_errors{4,1} = lm_errors;  
     %% Parametrization Errors and STD
    % 4 and 5GHz analysis
    if step <21
        pg_errors = struct('pg_err',zeros(1,3), 'pg_std',zeros(1,3), ...
            'pg_m_std',zeros(1,3),'pg_p_err',zeros(1,3),'pg_m_RSTD',...
            zeros(1,3));
        %mean error
        dif_pg_mv={n_p,1};                       
        m_dif_pg_mv=zeros(n_p,3);
        %percentual mean error 
        p_dif_pg_mv ={n_p,1};
        pm_dif_pg_mv =zeros(n_p,3);
         %standard deviation
        s_pg=zeros(n_p,3); RSTD_pg = zeros(1,3);
        for i=1:n_p
           %mean error
           dif_pg_mv{i,1} = abs(gen_A(i,1)-pg_A_nr{i,1});   
           m_dif_pg_mv(i,:) =  mean(dif_pg_mv{i,1});
           %percentual mean error
           p_dif_pg_mv{i,1} = (abs(gen_A(i,1)-pg_A_nr{i,1})...
               +gen_A(i,1))/gen_A(i,1)-1;
           pm_dif_pg_mv(i,:) = mean(p_dif_pg_mv{i,1});
           %standard deviation from true value of zero difference
           s_pg(i,:) = (sum(dif_pg_mv{i,1}).^2);
        end
        % Mean error
        mf = mean(m_dif_pg_mv);
        pg_errors.pg_err = mf;
    	% Percentual mean error
        pg_errors.pg_p_err= mean(pm_dif_pg_mv);
        % Standard deviation from zero 
        pg_errors.pg_std = sqrt(sum(s_pg)./(n_p*sg-1));
        for j=1:n_p
        %s_pg(j,1) = sum((m_dif_pg_mv(j,1)-dif_pg_mv{j,1}(:,1)).^2);
        %s_pg(j,2) = sum((m_dif_pg_mv(j,2)-dif_pg_mv{j,1}(:,2)).^2);
        %s_pg(j,3) = sum((m_dif_pg_mv(j,3)-dif_pg_mv{j,1}(:,3)).^2);    
        s_pg(j,1) = sum((mf(1,1)-dif_pg_mv{j,1}(:,1)).^2);
        s_pg(j,2) = sum((mf(1,2)-dif_pg_mv{j,1}(:,2)).^2);
        s_pg(j,3) = sum((mf(1,3)-dif_pg_mv{j,1}(:,3)).^2);        
        end
        std_pg=zeros(1,3);
        std_pg(1,:) = sqrt(sum(s_pg)/(n_p*sg-1));
        %RTSD
        RSTD_pg(1,1) = (std_pg(1,1)./mf(1,1))*100;
        RSTD_pg(1,2) = (std_pg(1,2)./mf(1,2))*100;
        RSTD_pg(1,3) = (std_pg(1,3)./mf(1,3))*100;
        % standard deviation form mean
        pg_errors.pg_m_std = std_pg;
        % Percentage relative standard deviation (RSTD) form mean
        pg_errors.pg_m_RSTD = RSTD_pg; 
        area_errors{3,1} = pg_errors;
    elseif step == 40 || step ==30
    % 2GHz calculation
        pg_errors = struct('pg_err',zeros(1,2), 'pg_std',zeros(1,2), ...
            'pg_m_std',zeros(1,2),'pg_p_err',zeros(1,2),'pg_m_RSTD',...
            zeros(1,2));
        %mean error
        dif_pg_mv={n_p,1};                       
        m_dif_pg_mv=zeros(n_p,2);
        %percentual mean error 
        p_dif_pg_mv ={n_p,1};
        pm_dif_pg_mv =zeros(n_p,2);
        %standard deviation
        s_pg=zeros(n_p,2); RSTD_pg = zeros(1,2);
        for i=1:n_p
           %mean error
           dif_pg_mv{i,1} = abs(gen_A(i,1)-pg_A2_nr{i,1});   
           m_dif_pg_mv(i,:) =  mean(dif_pg_mv{i,1});
           %percentual mean error
           p_dif_pg_mv{i,1} = (abs(gen_A(i,1)-pg_A2_nr{i,1})...
               +gen_A(i,1))/gen_A(i,1)-1;
           pm_dif_pg_mv(i,:) = mean(p_dif_pg_mv{i,1});
           %standard deviation from true value of zero difference
           s_pg(i,:) = (sum(dif_pg_mv{i,1}).^2);
        end
        % Mean error
        mf = mean(m_dif_pg_mv);
        pg_errors.pg_err = mf;
    	% Percentual mean error
        pg_errors.pg_p_err= mean(pm_dif_pg_mv);
        % Standard deviation from zero 
        pg_errors.pg_std = sqrt(sum(s_pg)./(n_p*sg-1));
        for j=1:n_p
        %s_pg(j,1) = sum((m_dif_pg_mv(j,1)-dif_pg_mv{j,1}(:,1)).^2);
        %s_pg(j,2) = sum((m_dif_pg_mv(j,2)-dif_pg_mv{j,1}(:,2)).^2);    
        s_pg(j,1) = sum((mf(1,1)-dif_pg_mv{j,1}(:,1)).^2);
        s_pg(j,2) = sum((mf(1,2)-dif_pg_mv{j,1}(:,2)).^2);
        end
        std_pg=zeros(1,2);
        std_pg(1,:) = sqrt(sum(s_pg)/(n_p*sg-1));
        %RTSD
        RSTD_pg(1,1) = (std_pg(1,1)./mf(1,1))*100;
        RSTD_pg(1,2) = (std_pg(1,2)./mf(1,2))*100;
        % standard deviation form mean
        pg_errors.pg_m_std = std_pg;
        % Percentage relative standard deviation (RSTD) form mean
        pg_errors.pg_m_RSTD = RSTD_pg;
        area_errors{3,1} = pg_errors;
    elseif step > 75
    %1GHz calculation    
        pg_errors = struct('pg_err',zeros(1,2), 'pg_std',zeros(1,2), ...
            'pg_m_std',zeros(1,2),'pg_p_err',zeros(1,2),'pg_m_RSTD',...
            zeros(1,2));
        %mean error
        dif_pg_mv={n_p,1};                       
        m_dif_pg_mv=zeros(n_p,1);
        %percentual mean errors
        p_dif_pg_mv ={n_p,1};
        pm_dif_pg_mv =zeros(n_p,1);
        %standard deviation
        s_pg=zeros(n_p,1); 
        for i=1:n_p
           %mean error
           dif_pg_mv{i,1} = abs(gen_A(i,1)-pg_A1_nr(i,:));   
           m_dif_pg_mv(i,:) =  mean(dif_pg_mv{i,1});
           %percentual mean error
           p_dif_pg_mv{i,1} = (abs(gen_A(i,1)-pg_A1_nr(i,:))...
               +gen_A(i,1))/gen_A(i,1)-1;
           pm_dif_pg_mv(i,:) = mean(p_dif_pg_mv{i,1});
           %standard deviation from true value of zero difference
           s_pg(i,:) = (sum(dif_pg_mv{i,1}).^2);
        end
        % Mean error
        mf = mean(m_dif_pg_mv);
        pg_errors.pg_err = mf;
    	% Percentual mean error
        pg_errors.pg_p_err= mean(pm_dif_pg_mv);
        % Standard deviation from zero 
        pg_errors.pg_std = sqrt(sum(s_pg)./(n_p*sg-1));   
        for j=1:n_p
            %s_pg(j,1) = sum((m_dif_pg_mv(j,1)-dif_pg_mv{j,1}(:,1)).^2);
            s_pg(j,1) = sum((mf-dif_pg_mv{j,1}(:,1)).^2);
        end
        std_pg = sqrt(sum(s_pg)/(n_p*sg-1));
        %RTSD
        RSTD_pg = (std_pg/mf(1,1))*100;
        % standard deviation form mean
        pg_errors.pg_m_std = std_pg;
        % Percentage relative standard deviation (RSTD) form mean
        pg_errors.pg_m_RSTD = RSTD_pg;
        area_errors{3,1} = pg_errors;
    end
end    