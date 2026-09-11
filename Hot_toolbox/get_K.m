function covariance_matrix = get_K(img_seq,avg)


covariance_matrix = zeros(length(avg), length(avg));
tmp_1 = zeros(1,length(img_seq));
tmp_2 = tmp_1;

for ii=1:1:length(img_seq)
    
    avg_seq{ii} = get_positive_values_only(img_seq{ii} - avg);
    
    avg_seq{ii} = avg_seq{ii}(90:170,90:170);
    
end

for xx=1:1:(length(avg_seq{1,1})^2)
    
    for ii=1:1:length(img_seq)
            
            tmp_1(ii) = avg_seq{ii}(xx);
    end
        
            for yy = 1:1:(length(avg_seq{1,1})^2)

                        for ii=1:1:length(img_seq)
            
                            tmp_2(ii) = avg_seq{ii}(yy);
    
                        end

                             covariance_matrix(xx,yy) = sum(tmp_1.*tmp_2) / length(tmp_1);
        
            end
    
    
    
end
