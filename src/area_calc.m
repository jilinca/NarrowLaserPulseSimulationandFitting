
function [areas_spl,areas_pol,areas_lm] = area_calc(all_mv_spl,all_mv_pol,...
    fwhm_spl,fwhm_pol,all_spl_y,all_spl_x,all_pol_y,all_pol_x,...
    all_mv_lm, all_fwhm_lm,n_p,sg)

 

 %% Method 1
 %areas_spl = all_mv_spl.*(fwhm_spl/(2*sqrt(2*log(2))))*sqrt(2*pi);
 %areas_pol = all_mv_pol.*(fwhm_pol/(2*sqrt(2*log(2))))*sqrt(2*pi);
 
 areas_lm = all_mv_lm.*(all_fwhm_lm/(2*sqrt(2*log(2))))*sqrt(2*pi);
 
 %% Method 2
 %%{
 areas_spl = zeros(n_p,sg); areas_pol = areas_spl;
 
 for i=1:n_p
     for j=1:sg
     areas_spl(i,:) = trapz(all_spl_x{i,1}(:,j), all_spl_y{i,1}(:,j));
     areas_pol(i,:) = trapz(all_pol_x{i,1}{j,1}, all_pol_y{i,1}{j,1}); 
     end
 end
 %%}
end