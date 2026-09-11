function fixed_pupil = fix_pupil_levels(pupil)

%Pupils from PAOLA are not ones and zeros but tend
%to taper off around the pupil edges. This short script
%fixes this

for xx=1:1:length(pupil)
       
       for yy=1:1:length(pupil)
            
            if (pupil(xx,yy) >=0.8)
     
                pupil(xx,yy) =1;
                
            elseif (pupil(xx,yy) <0.8)
                pupil(xx,yy) =0;
            end
       end
end

fixed_pupil=pupil;

end
       