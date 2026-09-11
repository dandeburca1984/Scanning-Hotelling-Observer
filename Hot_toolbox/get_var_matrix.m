%Calculate diagonal of K

function var_mat = get_var_matrix(cov_img)

matrix_size = size(cov_img);
N = matrix_size(1);
N_2 = matrix_size(2);
var_mat = zeros(1,length(cov_img));

for xx = 1:1:N_2
    
    tmp_1 = cov_img(:,xx);
    
    var_mat(xx) = sum(tmp_1.*tmp_1) / N;
    
end