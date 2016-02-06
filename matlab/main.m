%generate_MRF_Binomial([32,32],[-7.8 3.9 3.9])
generate_MRF([64,64],2);
%generate_MRF([32,32],5);
%generate_MRF_quick([32,32],2);
%generate_MRF_quick([32,32],5);
%ima = generate_MRF([32 32],2);
%denoise_MRF(ima);