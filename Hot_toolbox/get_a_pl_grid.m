function a_pl = get_a_pl_grid(r_pl,r_covariance,residual_H1,psf,star_position,A_star)
%Function to search for the optimal intensity of a companion

xshift = (star_position(1) - r_pl(1));
yshift = (star_position(2) - r_pl(2));
temp = shift_planet(psf,xshift,yshift,1);
r_temp = get_raster_image(temp);

sampling_1 = 1;
sampling_2 = 0.1;
index = 1;
a_pl_tmp1 = zeros(1,15);
f1 = a_pl_tmp1;
a_pl_tmp2 = a_pl_tmp1;
f2 = a_pl_tmp1;

for tmp_differential_mag = 1:sampling_1:15
    
    a_pl_tmp1(index) = A_star  / (10^(0.4*tmp_differential_mag));
    
    f1(index) = ( sum( ( (a_pl_tmp1(index)*r_temp) ./ r_covariance) ...
        .*(residual_H1 - 0.5*a_pl_tmp1(index)*r_temp)));
    
    index = index + 1;
    
end

[unused aa]=max(f1);
a_pl_ = a_pl_tmp1(aa);

delta_a_pl = 2.5*log10(A_star/a_pl_);

index=1;

for tmp_differential_mag2 = (delta_a_pl+sampling_1):-sampling_2:(delta_a_pl-sampling_1)
    
    a_pl_tmp2(index) = A_star  / (10^(0.4*tmp_differential_mag2));
    
    f2(index) = ( sum( ( (a_pl_tmp2(index)*r_temp) ./ r_covariance) ...
        .*(residual_H1 - 0.5*a_pl_tmp2(index)*r_temp)));
    
    index = index + 1;
    
end

[unused aa_]=max(f2);
a_pl = a_pl_tmp2(aa_);

end
    