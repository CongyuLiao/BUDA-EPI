function raw_corr = mrir_artifact_ghost_correct_CorrelationMethod(raw_roft, linear_fit_coeff,OS_factor)

% Kawin Setsompop
% 8/2/2012

% linear_fit_coeff here will be constant phase (row1) and secondlineshift
% (row2)

% input data needs to be in fourier space

if nargin == 2
    OS_factor = 5; % oversampling factor
end

if ( mrir_ice_dimensions(raw_roft, 'seg') < 2 ),
    error('uncorrected data contains only one segment');
end;


raw_roft = double(raw_roft);

% apply correction to forward data
fwd_lines = raw_roft(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:);
fwd_linesReshape = reshape(fwd_lines, size(fwd_lines,1), size(fwd_lines,2),[]);
fwd_linesReshapeCorrected = fwd_linesReshape;

raw_corr = raw_roft;

LineLength = size(raw_roft,1);
for DataSetCount = 1:size(fwd_linesReshape,3)
    Line1 = fwd_linesReshape(:,:,DataSetCount);
    ConstPhaseDiff = linear_fit_coeff(1,DataSetCount);
    SecondLineShiftBy = linear_fit_coeff(2,DataSetCount);
   
    PadOnEachSide = floor(LineLength*(OS_factor-1)/2); 
    Line1_OS = mrir_fDFT(padarray( mrir_iDFT(Line1,1),PadOnEachSide),1);
    Line1_OSshifted = circshift(Line1_OS,-SecondLineShiftBy);
    if(0)
        % use conjugate symetry approximation here.....
        if SecondLineShiftBy<0
            Line1_OSshifted(1:SecondLineShiftBy) = conj(Line1_OSshifted(1:SecondLineShiftBy));
        else
            Line1_OSshifted(end-SecondLineShiftBy:end) = conj(Line1_OSshifted(end-SecondLineShiftBy:end));
        end
    end
    
    Line1_OScorrected = Line1_OSshifted/exp(sqrt(-1)*ConstPhaseDiff);
    %fwd_linesReshapeCorrected(:,:,DataSetCount) = 3*Line1_OScorrected(2:3:end,:);  
    Image_Line1_OScorrected = mrir_iDFT(Line1_OScorrected,1);
    fwd_linesReshapeCorrected(:,:,DataSetCount) = mrir_fDFT(Image_Line1_OScorrected(1+PadOnEachSide:LineLength+PadOnEachSide,:),1);

end

% replace only the forward lines
dims = size(raw_corr(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:));
raw_corr(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:) = reshape(fwd_linesReshapeCorrected,dims);

raw_corrCollapsed = sum(raw_corr,8); 
raw_corrCollapsedUnCorrect = sum(raw_roft,8);

if (0)
    Slc = 48;
    for CoilCount =1:32
        figure(100); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsed(:,:,CoilCount,1,1,1,1,1,1,Slc))))));
        figure(101); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsedUnCorrect(:,:,CoilCount,1,1,1,1,1,1,Slc))))));
    end
    figure(102); imagesc(mrir_array_combine(mrir_image_crop(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsed(:,:,CoilCount,1,1,1,1,1,1,Slc))))),0));
    figure(103); imagesc(mrir_array_combine(mrir_image_crop(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsedUnCorrect(:,:,CoilCount,1,1,1,1,1,1,Slc))))),0));
    
end


  