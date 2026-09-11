function [test_stat position] = get_hot_grid(r_cov,r_data,PSF,star_location,a_pl,planet_positions)


test_stat_tmp = zeros(1,length(planet_positions));

for ii = 1:1:length(planet_positions)
    
     r_pl(1,1) = planet_positions(ii,1);
     r_pl(1,2) = planet_positions(ii,2);
     
     test_stat_tmp(ii) = hot_function_known_a_pl(r_pl,r_cov,r_data,PSF,star_location,a_pl);
         
end

[max_test_stat index] = max(test_stat_tmp);

test_stat = max_test_stat;

position(1,1) = planet_positions(index,1);
position(1,2) = planet_positions(index,2);

end