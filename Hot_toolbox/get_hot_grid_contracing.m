function [test_stat position] = get_hot_grid_contracing(r_cov,r_data,PSF,star_location,a_pl,planet_position_start)

sampling_1 = (1/2);
sampling_2 = (1/8);

index_x=1;
index = 1;

%First loop**************************
for rr= (planet_position_start(1,1)-1): sampling_1: (planet_position_start(1,1)+1)
    
    r_pl(1,1) = rr;
    index_y=1;
    
        for qq=(planet_position_start(1,2)-1) :sampling_1:(planet_position_start(1,2)+1)
            
            r_pl(1,2) = qq;
            
            positions(index,2) = qq;
            positions(index,1) = rr;
            
            test_stat_tmp(index_x,index_y) = hot_function_known_a_pl(r_pl,r_cov,r_data,PSF,star_location,a_pl);
            
            tmp_stat(index) = test_stat_tmp(index_x,index_y);
            
            index_y = index_y +1;
            
            index = index +1;
        end
        
        index_x=index_x + 1;
        
        
        
end

[max_of_image index ]=max(test_stat_tmp); 
[unused, location] = max(max_of_image);
test_stat_max = [location index(location)];
test_stat_max = test_stat_max-1;

position_tmp(1,1) =  (planet_position_start(1,1)-1) + ((test_stat_max(1,1))*sampling_1);
position_tmp(1,2) =  (planet_position_start(1,2)-1) + ((test_stat_max(1,2))*sampling_1);

[ee aa]=max(tmp_stat);
position_tmp_ = positions(aa,:);
%position_tmp(1,1) = planet_position_start(1,1) + ((test_stat_max(1,1) - round(length(test_stat_tmp)/2))*sampling_1  );
%position_tmp(1,2) = planet_position_start(1,2) + ((test_stat_max(1,2) - round(length(test_stat_tmp)/2))*sampling_1  );


%Second loop

index_x=1;
index=1;
for rr= (position_tmp(1,1)-(0.5) ) :sampling_2:(position_tmp(1,1)+(0.5))
    
    r_pl(1,1) = rr;
    index_y=1;
    
        for qq=(position_tmp(1,2)-(0.5) ):sampling_2:(position_tmp(1,2)+(0.5))
            
            r_pl(1,2) = qq;
            
            positions2(index,2) = qq;
            positions2(index,1) = rr;
            
            test_stat_tmp_prime(index_x,index_y) = hot_function_known_a_pl(r_pl,r_cov,r_data,PSF,star_location,a_pl);
            
            index_y = index_y +1;
            
            index = index +1;
        end
        
        index_x=index_x + 1;
        
end

[max_of_image index ]=max(test_stat_tmp_prime); 
[unused, location] = max(max_of_image);
position_ = [location index(location)] ;
%position_ = position_ - 1;

p_tmp = positions2( index(location),:);

position(1,1) = position_tmp(1,1)-(0.5) + (position_(1,1)*sampling_2);
position(1,2) = position_tmp(1,2)-(0.5) + (position_(1,2)*sampling_2);
%position(1,1) = position_tmp(1,1) + ((position_(1,1) - round(length(test_stat_tmp)/2))*sampling_2  );
%position(1,2) = position_tmp(1,2) + ((position_(1,2) - round(length(test_stat_tmp)/2))*sampling_2  );
%position_ = position_;
test_stat = test_stat_tmp_prime(position_(1,1),position_(1,2));


end
    