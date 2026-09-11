%% Test SCHO


clear;
clc;

data = fitsread('best_image_PSF38.fits');

psf = fitsread('best_image_PSF39.fits');

psf = make_gaussian(0,0,7,7,101);
psf2 = make_gaussian(15,-25,7,7,101);

data = psf*1e6 + psf2*1e5;

%results = Scanning_Hotelling_Observer(data,psf);

this  = Hoteling_Observer(data,psf);
this.Run_Observer


