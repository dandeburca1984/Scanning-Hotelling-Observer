function [ PSF ] = Phase2psf_beta( phase,pupil )

%phase = padarray(phase,[96 96],'pre');
%phase = padarray(phase,[95 95],'post');

%m = pupil;
%m = padarray(pupil,[96 96],'pre');
%m = padarray(m,[95 95],'post');
m=pupil;

amp=m.*complex(cos(phase),-sin(phase));

f=fft2(amp);
PSF=fftshift(abs(f).*abs(f));

end
