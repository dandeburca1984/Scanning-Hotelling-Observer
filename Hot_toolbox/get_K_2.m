function cov_mat = get_K_2(cov_img)


cov_mat= zeros(length(cov_img), length(cov_img));

matrix_size = size(cov_img);

N = matrix_size(1);
N_2 = matrix_size(2);

for xx = 1:1:N_2
    
    tmp_1 = cov_img(:,xx);
    
    for yy=1:1:N_2
        
        tmp_2 = cov_img(:,yy);
        
        cov_mat(xx,yy) = sum(tmp_1.*tmp_2) / N;
        
    end
    
end

