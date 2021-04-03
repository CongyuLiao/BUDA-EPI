function raw_corr = mrir_artifact_ghost_correct_CorrelationMethod_v2(raw_roft, linear_fit_coeff,OS_factor,ShiftFracForLine1)

% Kawin Setsompop
% 8/2/2012

% linear_fit_coeff here will be constant phase (row1) and secondlineshift
% (row2)

% input data needs to be in fourier space


if nargin == 2
    OS_factor = 5; % oversampling factor
    ShiftFracForLine1 = 0.5; % the fraction of the shift that will be applied to the odd line group
end

if ( mrir_ice_dimensions(raw_roft, 'seg') < 2 ),
    error('uncorrected data contains only one segment');
end;


raw_roft = double(raw_roft);

fwd_lines = raw_roft(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:);
fwd_linesReshape = reshape(fwd_lines, size(fwd_lines,1), size(fwd_lines,2),[]);
fwd_linesReshapeCorrected = fwd_linesReshape;

rev_lines = raw_roft(:,2:2:end,:,:,:,:,:,2,:,:,:,:,:,:,:,:);
rev_linesReshape = reshape(rev_lines, size(rev_lines,1), size(rev_lines,2),[]);
rev_linesReshapeCorrected = rev_linesReshape;

LineLength = size(raw_roft,1);
for DataSetCount = 1:size(fwd_linesReshape,3)
    Line1 = fwd_linesReshape(:,:,DataSetCount);
    Line2 = rev_linesReshape(:,:,DataSetCount);
    
    ConstPhaseDiff = linear_fit_coeff(1,DataSetCount);
    Total_SecondLineShiftBy = round(linear_fit_coeff(2,DataSetCount)); % can be non-interger due to median operation
    Line1Shift = floor(-Total_SecondLineShiftBy*ShiftFracForLine1);
    Line2Shift = Total_SecondLineShiftBy + Line1Shift;
   
    PadOnEachSide = floor(LineLength*(OS_factor-1)/2); 
    Line1_OS = mrir_fDFT(padarray( mrir_iDFT(Line1,1),PadOnEachSide),1);
    Line1_OSshifted = circshift(Line1_OS,Line1Shift);
    clear Line1_OS
    
    Line2_OS = mrir_fDFT(padarray( mrir_iDFT(Line2,1),PadOnEachSide),1);
    Line2_OSshifted = circshift(Line2_OS,Line2Shift);
    clear Line2_OS
    
    if(0)
        % use conjugate symetry approximation here.....
        if SecondLineShiftBy<0
            Line1_OSshifted(1:SecondLineShiftBy) = conj(Line1_OSshifted(1:SecondLineShiftBy));
        else
            Line1_OSshifted(end-SecondLineShiftBy:end) = conj(Line1_OSshifted(end-SecondLineShiftBy:end));
        end
    end
    
    Line1_OScorrected = Line1_OSshifted/exp(sqrt(-1)*ConstPhaseDiff);
 
    Image_Line1_OScorrected = mrir_iDFT(Line1_OScorrected,1);
    clear Line1_OScorrected
    Image_Line2_OScorrected = mrir_iDFT(Line2_OSshifted,1);
    clear Line1_shift
    
    fwd_linesReshapeCorrected(:,:,DataSetCount) = mrir_fDFT(Image_Line1_OScorrected(1+PadOnEachSide:LineLength+PadOnEachSide,:),1);
    clear Image_Line1_OScorrected
    rev_linesReshapeCorrected(:,:,DataSetCount) = mrir_fDFT(Image_Line2_OScorrected(1+PadOnEachSide:LineLength+PadOnEachSide,:),1);
    clear Image_Line2_OScorrected
end

raw_corr = raw_roft;
dims = size(raw_corr(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:));
raw_corr(:,1:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:) = reshape(fwd_linesReshapeCorrected,dims);
dims = size(raw_corr(:,2:2:end,:,:,:,:,:,1,:,:,:,:,:,:,:,:));
raw_corr(:,2:2:end,:,:,:,:,:,2,:,:,:,:,:,:,:,:) = reshape(rev_linesReshapeCorrected,dims);


if (0)
    raw_corrCollapsed = sum(raw_corr,8);
    raw_corrCollapsedUnCorrect = sum(raw_roft,8);
    Slc = 50;
    FigShift = 20; 
    
    % look at resulting images    
    figure(100+FigShift); 
    subplot(2,3,1); imagesc(mrir_array_combine(mrir_image_crop(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsedUnCorrect(:,:,:,1,1,1,1,1,1,Slc))))),0)); title('Uncorrected')
    subplot(2,3,2); imagesc(mrir_array_combine(mrir_image_crop(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsed(:,:,:,1,1,1,1,1,1,Slc))))),0)); title('Corrected')
    subplot(2,3,3); imagesc(mrir_array_combine(mrir_image_crop(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsed(:,1:2:end,:,1,1,1,1,1,1,Slc))))),0)); title('Odd line')
  
    for CoilCount =1:32
        figure(101+FigShift); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsedUnCorrect(:,:,CoilCount,1,1,1,1,1,1,Slc))))),[0 10]); title('Uncorrected')
        figure(102+FigShift); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsed(:,:,CoilCount,1,1,1,1,1,1,Slc))))),[0 10]); title('Corrected')
        figure(103+FigShift); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsedUnCorrect(:,1:2:end,CoilCount,1,1,1,1,1,1,Slc))))),[0 10]); title('Oddlines');
        %figure(104); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsed(:,2:2:end,CoilCount,1,1,1,1,1,1,Slc))))),[0 50]);        
        %figure(104); subplot(5,7,CoilCount); imagesc(abs(mrir_iDFT_phasencode(mrir_iDFT_freqencode(double(raw_corrCollapsedUnCorrect(:,:,CoilCount,1,1,1,1,1,1,Slc))))),[0 50]);

    end
    
    % look at k-space data
    
    CoilToLookAt = [7, 9, 23,27]; 
    centerKxline = size(raw_corrCollapsed,1) *(1/2)   +1; 
    centerKyline = size(raw_corrCollapsed,2) *(1/3)   +1; % 6/8 PF
    windowKx = -ceil(size(raw_corrCollapsed,1) *(1/8)):ceil(size(raw_corrCollapsed,1) *(1/8)) ;
    windowKy = -ceil(size(raw_corrCollapsed,2) *(1/8)):ceil(size(raw_corrCollapsed,2) *(1/8)) ;
    
    for CoilCount = 1:length(CoilToLookAt)
        figure(201+FigShift); subplot(2,3,CoilCount); imagesc((abs(raw_corrCollapsed(windowKx+centerKxline,windowKy+centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc))));
        %figure(202); subplot(2,3,CoilCount); imagesc((abs(raw_corrCollapsedUnCorrect(windowKx+centerKxline,windowKy+centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc))));
        
        figure(203+FigShift); subplot(2,3,CoilCount); imagesc(angle(raw_corrCollapsed(windowKx+centerKxline,windowKy+centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));
        figure(204+FigShift); subplot(2,3,CoilCount); imagesc(angle(raw_corrCollapsedUnCorrect(windowKx+centerKxline,windowKy+centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));
                
%         figure(205); subplot(2,3,CoilCount); plot(abs(raw_corrCollapsed(windowKx+centerKxline,centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));
%         hold on; subplot(2,3,CoilCount); plot(abs(raw_corrCollapsed(windowKx+centerKxline,centerKyline+1,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)),'r'); 
%         
%         figure(206); subplot(2,3,CoilCount); plot(angle(raw_corrCollapsed(windowKx+centerKxline,centerKyline,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));
%         hold on; subplot(2,3,CoilCount); plot(angle(raw_corrCollapsed(windowKx+centerKxline,centerKyline+1,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)),'r'); 
%                 
%          figure(207); subplot(2,3,CoilCount); plot(abs(raw_corrCollapsed(centerKxline,:,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));       
%          figure(208); subplot(2,3,CoilCount); plot(angle(raw_corrCollapsed(centerKxline,:,CoilToLookAt(CoilCount),1,1,1,1,1,1,Slc)));

    end
     
end


  