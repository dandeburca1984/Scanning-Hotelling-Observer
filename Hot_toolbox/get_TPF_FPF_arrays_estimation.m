function [FPF_array, TPF_array] = get_TPF_FPF_arrays_estimation(test_values,tol)

num_true_locations = size(test_values, 1);

N = size(test_values, 2);

values_planet_present = zeros(num_true_locations * N, 1);

values_planet_absent = zeros(num_true_locations * N, 1);

values_distance = zeros(num_true_locations * N, 1);

values_differential_magnitude = zeros(num_true_locations * N, 1);

FPF_array = zeros(1,N*2);
TPF_array = zeros(1,N*2);

index = 1;


    for ii = 1:N
        
        values_planet_present(index) = test_values{ii}.value_planet_present;
        
        values_planet_absent(index) = test_values{ii}.value_planet_absent;
        
        values_distance(index) = test_values{ii}.distance;
        
        values_differential_magnitude(index) = test_values{ii}.mag_error; 
        
        index = index + 1;
        
    end


%disp('values_distance');

sorted_values = sort([values_planet_present; values_planet_absent]);

sorted_values = [sorted_values(1) - 1; sorted_values; sorted_values(end) + 1];


for ii = 1:(2 * num_true_locations * N )
    
    %t = sorted_values(ii);
    t = sorted_values(2 * N * num_true_locations + 1 - ii);
    
    FPF_array(ii) = sum(values_planet_absent > t) / (N*num_true_locations);
    
    TPF_array(ii) = sum((values_planet_present >= t) & (values_distance <= tol(1,1)) & (values_differential_magnitude<=tol(1,2)))   / (N*num_true_locations);
    
     
end



end
    

