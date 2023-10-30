 function [out] = ELUM(img)
% ELUM (energy-weighted local uncertainty measure)
% 
% 1.Please note that this code do not contain the segmentation step. If you
% want to get the locations of small targets, you need to use some
% segmentation algorithm.
% 2.Please ensure that the input is a single channel and that the data type
% is double. 
%==========================================================================
% If you use this code in your publications, please cite:
% E. Zhao, W. Zheng, M. Li, H. Sun and J. Wang, "Infrared Small Target
% Detection Using Local Component Uncertainty Measure With Consistency
% Assessment," in IEEE Geoscience and Remote Sensing Letters, vol. 19, pp.
% 1-5, 2022, Art no. 6518205, doi: 10.1109/LGRS.2022.3221088.   
%==========================================================================

    batchsize = 3;

    mean_krl = ones(batchsize, batchsize) / batchsize^2;
    I_mean = imfilter(img,mean_krl,'replicate'); 
    max_krl = ones(batchsize,batchsize);
    I_max = ordfilt2(img,batchsize^2,max_krl,'symmetric'); 

    I_LC =( img ./ I_mean ) .* ( img ./ I_max) ; % Component Consistency Evaluation

    % gauss_krl = fspecial('gaussian',[ batchsize , batchsize ],1); 
    gauss_krl=[1,2,1;2,4,2;1,2,1]./16;
    I_gauss=imfilter(img,gauss_krl,'replicate');

    %% Energy Weighting Factor
    backMean_krl = ones(batchsize,batchsize); 
    backMean_krl( round( batchsize/2), round( batchsize/2)) = 0;
    backMean_krl = backMean_krl/(batchsize^2-1);
    I_backMean = imfilter(I_gauss,backMean_krl,'replicate');
    I_Gain = max ( img - I_backMean , 0 ) ;  
    %% Local Uncertainty Measurement
    minEntory = - 9 * (1/9)*log((1-1/9).^2);         % minimum entropy
    LcBlock = I_LC ./ imfilter(I_LC,ones(3,3),'replicate'); 
    Entory = LcBlock .* log( (1-LcBlock).^2 ); 
    Entory(LcBlock<=0) = 0;
    I_Entory = imfilter(Entory,-ones(3,3),'replicate') - minEntory ;  % Uncertainty Measurement

 
    I_WEntory = I_Gain  .* I_Entory;
    I_WEntory = max(I_WEntory , 0 ); 
    out = ( I_WEntory - min(I_WEntory(:)) ) ./ (max(I_WEntory(:)) - min(I_WEntory(:))) .*255;    
end
