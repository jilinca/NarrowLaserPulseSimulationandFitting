function [fwhm_spl, fwhm_pol]=...
    fwhm_calc(all_spl_y,all_spl_x,all_pol_y,all_pol_x,n_p,sg)
fwhm_spl = zeros(n_p,sg); fwhm_pol= fwhm_spl;
%% Calculation for spline and polynomial
for i=1:n_p
    for j=1:sg
        mv = max(all_spl_y{i,1}(:,j));
        if length(mv)>1
            mv = mv(1);
        end
        %% Spline
        hm = mv/2;
        idx1 = find(all_spl_y{i,1}(:,j)>hm,1) +[-1 0];
        idx2 = find(all_spl_y{i,1}(:,j)>hm,1,'last') + [0 1];
        x1 = interp1(all_spl_y{i,1}(idx1,j),all_spl_x{i,1}(idx1,j),hm);
        x2 = interp1(all_spl_y{i,1}(idx2,j),all_spl_x{i,1}(idx2,j),hm);
        fwhm_spl(i,j)= x2-x1;
        %% Polynomial
        mv_pol = max(all_pol_y{i,1}{j,1}(1,:));
        if length(mv_pol)>1
            mv_pol=mv_pol(1);
        end
        hm_pol = mv_pol/2;
        idx1 = find(all_pol_y{i,1}{j,1}(1,:)>hm_pol,1)+[-1 0];
        % if clause for 1GHz case
        %%{
        if idx1(1,1) ==0
            idx1(1,1) =1;
            idx1(1,2) = 2;
        end
        %%}
        idx2 = find(all_pol_y{i,1}{j,1}(1,:)>hm_pol,1,'last')+[0 1];
        % if clause for 1GHz case        
        %%{
        if length(all_pol_y{i,1}{j,1}) == idx2(1,1)
           idx2(1,2) = idx2(1,1);
           idx2(1,1) = idx2(1,1)-1;
        end
        x1 = interp1(all_pol_y{i,1}{j,1}(1,idx1),...
            all_pol_x{i,1}{j,1}(1,idx1),hm_pol);
        x2 = interp1(all_pol_y{i,1}{j,1}(1,idx2),...
            all_pol_x{i,1}{j,1}(1,idx2),hm_pol);
        fwhm_pol(i,j) = abs(x2-x1);
    end
end    
end