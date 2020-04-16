% Calculates error for noise removed functions
function [peak_value_errors,difs_for_hists,p_nums] =peak_err(...
    peak_value_errors, gen_peaks,all_mv_smp,all_mv_pol,all_mv_spl,...
    all_mv_lm,pg_mv_nr,pg_mv2_nr,pg_mv1_nr, n_p,sg,step)
    %% Cell for storing histogram data
    
    if step <21
        difs_for_hists ={7,1};
    elseif step == 30 || step == 40
        difs_for_hists ={6,1};
    else
        difs_for_hists = {5,1};
    end    
        
    %Temporary matrix for storing histogram data
    dif_h = zeros(n_p,sg); dif_h1 =dif_h; dif_h2 =dif_h; dif_h3 =dif_h;
    %% RETREIVED PEAK VALUE ERRORS AND STD
    smp_errors = struct('smp_err',0, 'smp_std',0, 'smp_m_std',0,...
    'smp_p_err',0,'smp_m_RSTD',0);
    
    %mean error
    dif_smp_mv={n_p,1};                       
    m_dif_smp_mv =zeros(n_p,1);
    %percentage mean error
    p_dif_smp_mv = {n_p,1};                   
    pm_dif_smp_mv =zeros(n_p,1);
    %standard deviation
    s_smp=zeros(n_p,1);
    
    for i=1:n_p
       %mean error
       dif_smp_mv{i,1} = abs(gen_peaks(i,1)-all_mv_smp(i,:));    
       dif_h(i,:) = gen_peaks(i,1)-all_mv_smp(i,:);             
       m_dif_smp_mv(i,1) = mean(dif_smp_mv{i,1}(1,:));                    
       %percentage mean error
       p_dif_smp_mv{i,1} = (abs(gen_peaks(i,1)-all_mv_smp(i,:))...
           +gen_peaks(i,1))/gen_peaks(i,1) -1;
       pm_dif_smp_mv(i,1) = mean(p_dif_smp_mv{i,1}(1,:));
       % Standard deviation from true value of zero differnece
       s_smp(i,1) = sum((dif_smp_mv{i,1}(1,:)).^2);
    end
    
    difs_for_hists{1,1} = reshape(dif_h,[n_p*sg,1]); 
    
    %Mean Error
    mf_smp =  mean(m_dif_smp_mv(:,1));  
    smp_errors.smp_err = mf_smp;                                           
    % Percentual mean error
    smp_errors.smp_p_err = mean(pm_dif_smp_mv(:,1)); 
    % Standard deviation from true value of zero difference
    smp_errors.smp_std = sqrt(sum(s_smp)/(sg*n_p -1));
    % Standard deviation from mean differences
    for j=1:n_p
        %s_smp(i,1) = sum((m_dif_smp_mv(i,1)-dif_smp_mv{i,1}(1,:)).^2);
        s_smp(i,1) = sum((mf_smp-dif_smp_mv{i,1}(1,:)).^2);
    end
    std_smp =sqrt(sum(s_smp)/(sg*n_p-1));
    %RTSD
    RSTD_smp = (std_smp/mf_smp)*100;
    
    smp_errors.smp_m_std = sqrt(sum(s_smp)/(sg*n_p -1));
    % Percentage or relative standard deviation (RSTD)
    smp_errors.smp_m_RSTD = RSTD_smp;
    
    peak_value_errors{1,1} = smp_errors;
    
   
    %% SPLINE ERRORS AND STD
    spl_errors = struct('spl_err',0, 'spl_std',0, 'spl_m_std',0,...
        'spl_p_err',0,'spl_m_RSTD',0);
    
    %mean error
    dif_spl_mv={n_p,1};                       
    m_dif_spl_mv =zeros(n_p,1);
    %percentage mean error
    p_dif_spl_mv = {n_p,1};                   
    pm_dif_spl_mv =zeros(n_p,1);
    %standard deviation
    s_spl=zeros(n_p,1); 
    
    for i=1:n_p
       %mean error
       dif_spl_mv{i,1} = abs(gen_peaks(i,1)-all_mv_spl(i,:));
       dif_h(i,:) = gen_peaks(i,1)-all_mv_spl(i,:);
       m_dif_spl_mv(i,1) = mean(dif_spl_mv{i,1}(1,:));                      
       %percentage mean error
       p_dif_spl_mv{i,1} = (abs(gen_peaks(i,1)-all_mv_spl(i,:))...
           +gen_peaks(i,1))/gen_peaks(i,1) -1;
       pm_dif_spl_mv(i,1) = mean(p_dif_spl_mv{i,1}(1,:));
       % Standard deviation from true value of zero differnece
       s_spl(i,1) = sum((dif_spl_mv{i,1}(1,:)).^2);
    end
    
    difs_for_hists{2,1} = reshape(dif_h,[n_p*sg,1]); 
    
    %Mean Error
    spl_mf =  mean(m_dif_spl_mv(:,1));  
    spl_errors.spl_err = spl_mf;                                            
    % Percentual mean error
    spl_errors.spl_p_err = mean(pm_dif_spl_mv(:,1)); 
    % Standard deviation from true value of zero difference
    spl_errors.spl_std = sqrt(sum(s_spl)/(sg*n_p -1));
    % Standard deviation from mean differences
    for j=1:n_p
        %s_spl(j,1) = sum((m_dif_spl_mv(j,1)-dif_spl_mv{j,1}(1,:)).^2);
        s_spl(j,1) = sum((spl_mf-dif_spl_mv{j,1}(1,:)).^2);
    end
    
    std_spl = sqrt(sum(s_spl)/(sg*n_p-1));
    RSTD_spl= (std_spl/spl_mf)*100;
    
    spl_errors.spl_m_std = std_spl;
    % Percentage or relative standard deviation (RSTD)
    spl_errors.spl_m_RSTD = RSTD_spl;
    
    peak_value_errors{2,1} = spl_errors;
    
    %% POLYNOMIAL ERRORS AND STD
    % KEY: 
    % poly_err describes the mean error from generated values
    % poly std describes the standard deviation from a mean difference of
    % zero
    % poly_m_std describs the standard deviation from the mean error of
    % generated values
    % poly_p_err displays mean error from generated values as a percentage
    % poly_m_RSTD describes the standard deviation from the mean error 
    %of generated values as a percentage
    % poly_p_std describes the standard deviation from a mean difference
    % error of zero
    
    poly_errors = struct('pol_err',0, 'pol_std',0, 'pol_m_std',0,...
        'pol_p_err',0,'pol_m_RSTD',0);
    %mean error
    dif_pol_mv={n_p,1};                       
    m_dif_pol_mv =zeros(n_p,1);
    %percentage mean error
    p_dif_pol_mv = {n_p,1};                   
    pm_dif_pol_mv =zeros(n_p,1);
    %standard deviation
    s_pol=zeros(n_p,1);
    
    for i=1:n_p
        %mean error
        dif_pol_mv{i,1} = abs(gen_peaks(i,1)-all_mv_pol(i,:));    
        dif_h(i,:) = gen_peaks(i,1)-all_mv_pol(i,:);
        m_dif_pol_mv(i,1) = mean(dif_pol_mv{i,1}); 
        %percentage mean error
        p_dif_pol_mv{i,1} = (abs(gen_peaks(i,1)-all_mv_pol(i,:))...
            +gen_peaks(i,1))/gen_peaks(i,1) -1;
        pm_dif_pol_mv(i,1) = mean(p_dif_pol_mv{i,1}(1,:));
        % Standard deviation from true value of zero differnece
        s_pol(i,1) = sum((dif_pol_mv{i,1}(1,:)).^2);
    end
    
    difs_for_hists{3,1} = reshape(dif_h,[n_p*sg,1]);
    
    %Mean Error
    mf_pol =  mean(m_dif_pol_mv(:,1));  
    poly_errors.pol_err = mf_pol;     % mean difference from zero 
    % Percentual mean error
    poly_errors.pol_p_err = mean(pm_dif_pol_mv(:,1)); 
    % Standard deviation from true value of zero difference
    poly_errors.pol_std = sqrt(sum(s_pol)/(sg*n_p -1));
    % Standard deviation from mean differences
    
    for j=1:n_p
        %s_pol(j,1) = sum((m_dif_pol_mv(j,1)-dif_pol_mv{j,1}(1,:)).^2);
        s_pol(j,1) = sum((mf_pol-dif_pol_mv{j,1}(1,:)).^2);
    end
    
    std_pol = sqrt(sum(s_pol(:,1))/(sg*n_p-1));
    RSTD_pol= (std_pol/mf_pol)*100;
    
    poly_errors.pol_m_std = std_pol;
    % Percentage or relative standard deviation (RSTD)
    poly_errors.pol_m_RSTD = RSTD_pol;
    
    peak_value_errors{3,1} = poly_errors;
    
    %% LEVENBERG MARQUARDT ERRORS 
    lm_errors = struct('lm_err',0, 'lm_std',0, 'lm_m_std',0,...
        'lm_p_err',0,'lm_m_RSTD',0);
    %mean error
    dif_lm_mv=zeros(1,sg);                       
    m_dif_lm_mv =zeros(1,1);
    %percentage mean error
    p_dif_lm_mv = zeros(1,sg);                   
    pm_dif_lm_mv =zeros(1,1);
    %standard deviation
    s_lm=zeros(1,1);
    %for histogram
    dif_hlm = zeros(1,sg);
    
    p_nums =[];
    
    ANoP =0;
    % Not all the pulses are anaylzed in case of the LM algorithm since 
    % sometimes the algorithm doesn not find the good parmaeters and just 
    % produces the intial guess again
    for i=1:n_p
        if sum(all_mv_lm(i,:)) ~= 5*600
            if ANoP == 0
                % mean error
                dif_lm_mv(i,:) = abs(gen_peaks(i,1)-all_mv_lm(i,:));
                dif_hlm(i,:) = gen_peaks(i,1)-all_mv_lm(i,:);
                m_dif_lm_mv(i,1) = mean(dif_lm_mv(i,:));
                %percentage mean error
                p_dif_lm_mv(i,:) = (abs(gen_peaks(i,1)-all_mv_lm(i,:))...
                    +gen_peaks(i,1))/gen_peaks(i,1) -1;
                pm_dif_lm_mv(i,1) = mean(p_dif_lm_mv(i,:));
                %preliminary sd calc from true value of zero
                s_lm(i,1) = sum((dif_lm_mv(1,:)).^2);
                ANoP = ANoP +1; 
                p_nums=[p_nums, i];
            end    
            if ANoP >=1
                ANoP = ANoP +1;                                              % Number of pulses analyzed
                %mean error
                dif_lm_mv = [dif_lm_mv;abs(gen_peaks(i,1)-all_mv_lm(i,:))]; %#ok<AGROW>
               
                dif_hlm = [dif_hlm; gen_peaks(i,1)-all_mv_lm(i,:)];         %#ok<AGROW>
                m_dif_lm_mv =[m_dif_lm_mv; mean(dif_lm_mv(ANoP,:))];        %#ok<AGROW>
                %percentage mean error
                p_dif_lm_mv = [p_dif_lm_mv;...                              
                    ((abs(gen_peaks(i,1)-all_mv_lm(i,:))...
                    +gen_peaks(i,1))/gen_peaks(i,1) -1)];                   %#ok<AGROW>
                pm_dif_lm_mv = [pm_dif_lm_mv; mean(p_dif_lm_mv(ANoP,:))];   %#ok<AGROW>
                %preliminary sd calc from true value of zero
                s_lm = [s_lm; sum((dif_lm_mv(ANoP,:)).^2)];                 %#ok<AGROW>
                p_nums=[p_nums,i];
            end
        end
    end
    
    
    difs_for_hists{4,1} = reshape(dif_hlm,[length(dif_hlm(:,1))*sg,1]);
    
    %Mean Error
    mf_lm =  mean(m_dif_lm_mv(:,1));  
    lm_errors.lm_err = mf_lm;                  % mean difference from zero error 
    % Percentual mean error
    lm_errors.lm_p_err = mean(pm_dif_lm_mv(:,1)); 
    % Standard deviation from true value of zero difference
    lm_errors.lm_std = sqrt(sum(s_lm)/(sg*ANoP -1));
    % Standard deviation from mean differences
    s_lm = zeros(ANoP,1);
    for j=1:ANoP
        %s_lm(j,1) = sum((m_dif_lm_mv(j,1)-dif_lm_mv{j,1}(1,:)).^2);
        s_lm(j,1) = sum((mf_lm-dif_lm_mv(j,:)).^2);
    end
    
    std_lm = sqrt(sum(s_lm(:,1))/(sg*ANoP-1));
    RSTD_lm= (std_lm/mf_lm)*100;
    
    lm_errors.lm_m_std = std_lm;
    % Percentage or relative standard deviation (RSTD)
    lm_errors.lm_m_RSTD = RSTD_lm;
    
    peak_value_errors{5,1} = lm_errors;
    
    %% PARAMETRIZATION ERRORS and STD
    
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
           dif_pg_mv{i,1} = abs(gen_peaks(i,1)-pg_mv_nr{i,1});
           % For histogram
           dif_h1(i,:) = gen_peaks(i,1)-pg_mv_nr{i,1}(:,1);
           dif_h2(i,:) = gen_peaks(i,1)-pg_mv_nr{i,1}(:,2);
           dif_h3(i,:) = gen_peaks(i,1)-pg_mv_nr{i,1}(:,3);
           
           m_dif_pg_mv(i,:) =  mean(dif_pg_mv{i,1});
           %percentual mean error
           p_dif_pg_mv{i,1} = (abs(gen_peaks(i,1)-pg_mv_nr{i,1})...
               +gen_peaks(i,1))/gen_peaks(i,1)-1;
           pm_dif_pg_mv(i,:) = mean(p_dif_pg_mv{i,1});
           %standard deviation from true value of zero difference
           s_pg(i,:) = (sum(dif_pg_mv{i,1}).^2);
        end
       
        %For histogram
        difs_for_hists{5,1} = reshape(dif_h1,[n_p*sg,1]);
        difs_for_hists{6,1} = reshape(dif_h2,[n_p*sg,1]);
        difs_for_hists{7,1} = reshape(dif_h3,[n_p*sg,1]);
        
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
        
        peak_value_errors{4,1} = pg_errors;
    
    elseif step == 40 || step == 30
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
           dif_pg_mv{i,1} = abs(gen_peaks(i,1)-pg_mv2_nr{i,1});   
           m_dif_pg_mv(i,:) =  mean(dif_pg_mv{i,1});
           
           %For Histograms
           dif_h1(i,:) = gen_peaks(i,1)-pg_mv_nr{i,1}(:,1);
           dif_h2(i,:) = gen_peaks(i,1)-pg_mv_nr{i,1}(:,2);
           
           %percentual mean error
           p_dif_pg_mv{i,1} = (abs(gen_peaks(i,1)-pg_mv2_nr{i,1})...
               +gen_peaks(i,1))/gen_peaks(i,1)-1;
           pm_dif_pg_mv(i,:) = mean(p_dif_pg_mv{i,1});
           %standard deviation from true value of zero difference
           s_pg(i,:) = (sum(dif_pg_mv{i,1}).^2);
        end
        
        %For Histogram plotting 
        difs_for_hists{5,1} = reshape(dif_h1,[n_p*sg,1]);
        difs_for_hists{6,1} = reshape(dif_h2,[n_p*sg,1]);
       
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
        
        peak_value_errors{4,1} = pg_errors;
        
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
           dif_pg_mv{i,1} = abs(gen_peaks(i,1)-pg_mv1_nr(i,:)); 
           dif_h1(i,:) = gen_peaks(i,1)-pg_mv_nr{i,1}(:,1);
           m_dif_pg_mv(i,:) =  mean(dif_pg_mv{i,1});
           %percentual mean error
           p_dif_pg_mv{i,1} = (abs(gen_peaks(i,1)-pg_mv1_nr(i,:))...
               +gen_peaks(i,1))/gen_peaks(i,1)-1;
           pm_dif_pg_mv(i,:) = mean(p_dif_pg_mv{i,1});
           %standard deviation from true value of zero difference
           s_pg(i,:) = (sum(dif_pg_mv{i,1}).^2);
        end
        
        %For Histogram plotting 
        difs_for_hists{4,1} = reshape(dif_h1,[n_p*sg,1]);
       
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
        
        peak_value_errors{4,1} = pg_errors;
    end
end