function [FPF_array, TPF_array] = get_TPF_FPF_arrays(test_values,tol)

num_true_locations = size(test_values, 1);

N = size(test_values, 2);

values_planet_present = zeros(num_true_locations * N, 1);

values_planet_absent = zeros(num_true_locations * N, 1);

values_distance = zeros(num_true_locations * N, 1);

index = 1;


    
    for ii = 1:N
        
        values_planet_present(index) = test_values{ii}.value_planet_present;
        
        values_planet_absent(index) = test_values{ii}.value_planet_absent;
        
        values_distance(index) = test_values{ii}.distance;
        
        index = index + 1;
        
    end


%disp('values_distance');

sorted_values = sort([values_planet_present; values_planet_absent]);

sorted_values = [sorted_values(1) - 1; sorted_values; sorted_values(end) + 1];

%disp('sorted_values');

FPF_array = [0];
TPF_array = [0];
FP = -1;
FN = -1;
TP=1;

for ii = 1:(2 * num_true_locations * N )
    
    %t = sorted_values(ii);
    t = sorted_values(2 * N * num_true_locations + 1 - ii);
    
    tmp_FP = sum(values_planet_absent > t);
    
    tmp_TP = sum((values_planet_present >= t) & (values_distance <= tol));
    
    %tmp_FN = sum(values_planet_present <= t);
    
    
    %if (tmp_FP ~= FP) &&(tmp_FN ~= FN)%(tmp_TP~=TP)
    if  (tmp_FP ~= FP) &&(tmp_TP~=TP)
    
        FP = tmp_FP;
       % FN = tmp_FN;
        
        TP = tmp_TP;
        %TP = (N*num_true_locations) - FN;
        %TN = (N*num_true_locations) - FP;
        
        % Update the arrays
        TPF_array = [TPF_array, TP / (N*num_true_locations)];
        FPF_array = [FPF_array, FP / (N*num_true_locations)];
        
        %FPF_array = [FP / (N * num_true_locations), FPF_array];
        
        %TPF_array = [TP / (N * num_true_locations), TPF_array];
        
    end
end

% Make sure that the plot ends at (1, 1)
%TPF_array = [TPF_array, 1];
%FPF_array = [FPF_array, 1];

end

