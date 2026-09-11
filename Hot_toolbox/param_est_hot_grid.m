function [test_stat position int] = param_est_hot_grid(r_cov,r_data,PSF,star_location,planet_position_start,A_star)

sampling_1 = 1;
sampling_2 = (1/20);

index = 1;
grid_size_1 = 3;
start_location =  star_location - planet_position_start; 

%First loop**************************
for rr= (start_location(1,1)-grid_size_1): sampling_1: (start_location(1,1)+grid_size_1)
   
    r_pl(1,1) = rr;
   
        for qq=(start_location(1,2)-grid_size_1) :sampling_1:(start_location(1,2)+grid_size_1)
            
            r_pl(1,2) = qq;
            
            positions(index,1) = rr;
            positions(index,2) = qq;
            
            a_pl_tmp(index) = get_a_pl_grid(r_pl,r_cov,r_data,PSF,star_location,A_star);
            
            tmp_stat(index) = hot_function_known_a_pl(r_pl,r_cov,r_data,PSF,star_location,a_pl_tmp(index));
            
            index = index +1;
        end
        
   
end

[unused aa]=max(tmp_stat);
position_ = positions(aa,:);


index = 1;

%First loop**************************
for rr= (position_(1,1)-sampling_1): sampling_2: (position_(1,1)+sampling_1)
    
    r_pl(1,1) = rr;
 
        for qq=(position_(1,2)-sampling_1) :sampling_2:(position_(1,2)+sampling_1)
            
            r_pl(1,2) = qq;
            
            positions_(index,1) = rr;
            positions_(index,2) = qq;
            
            a_pl_tmp(index) = get_a_pl_grid(r_pl,r_cov,r_data,PSF,star_location,A_star);
           
            tmp_stat_(index) = hot_function_known_a_pl(r_pl,r_cov,r_data,PSF,star_location,a_pl_tmp(index));
            
            index = index +1;
        end
       
end

[test_stat aa]=max(tmp_stat_);
position = positions_(aa,:) -  star_location;
int = a_pl_tmp(aa);


end