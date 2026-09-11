function retval = get_AUC(FPF_array, TPF_array)
retval = abs(sum((FPF_array(2:end) - FPF_array(1:(end - 1))) .* (TPF_array(2:end) + TPF_array(1:(end - 1))))) / 2;

