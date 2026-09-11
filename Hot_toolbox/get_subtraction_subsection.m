%function to calculate the subtraction subsection

function Sr = get_subtraction_subsection(ck,K,rmin,theta_min,theta_max,image_centre,dr)

Sr = zeros(size(K{1}));
count=0;

for qq=1:1:length(K)
    
    rmax = rmin+dr;
    
    [Sk ind] = get_image_slice_and_pixels(K{qq},rmin,rmax,theta_min,theta_max,image_centre);
    
    tmp = ck(qq) .* Sk;
    
    count = count + Sk;
    
    Sr = Sr +tmp;
    
end
