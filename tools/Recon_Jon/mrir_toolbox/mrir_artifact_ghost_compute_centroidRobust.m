function linear_fit_coeff = mrir_artifact_ghost_compute_centroidRobust(phascor,PercentDataUse,RegressionMethod)
%MRIR_ARTIFACT_GHOST_COMPUTE
%
% linear_fit_coeff = mrir_artifact_ghost_compute(phascor_raw)
%
% returns a 2 x Nphascor matrix of linear fit coefficients, where Nphascor
% is the number of unique phase correction lines collected at the center of
% k-space (e.g., the number of coil channels). at least one pair of phase
% correction navigators must be collected, and hopefully this measurement is
% repeated several times. the first row is the offset term, and the second
% row is the slope. all units are radians.

% NOTE: for single-shot (unsegmented) EPI, the number of segments of the
% phascor data (i.e., NSeg) should be 1, for compatibility. unfortunately,
% the convention for the raw data is that the two segments labeled by the
% loopcounters remain, so, for the raw data only, NSeg = 2.

% NOTE: algorithm is based on fitting a 1st order polynomial to the phase
% difference between two navigators in hybrid space where the zeroth order
% term represents the global phase changes attributable to global B0
% differences between the two navigators, and the first order term
% represents the spatially varying phase changes in hybrid space which
% correspond to a shift in the echo in k space. the typical method for
% calculating these two terms is based on an autocorrelation of the
% navigators, introduced by Ahn & Cho (1987) and summarized by Clare in
% equation (4.3),
% http://users.fmrib.ox.ac.uk/~stuart/thesis/chapter_4/section4_3.html

% references:
%
%  Bruder H, Fischer H, Reinfelder HE, Schmitt F. Image reconstruction for
%  echo planar imaging with nonequidistant k-space sampling. Magn Reson Med.
%  1992 Feb;23(2):311-23. PMID: 1549045.


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2006/dec/31
% $Id: mrir_artifact_ghost_compute.m,v 1.6 2008/05/16 17:58:36 jonnyreb Exp $
%**************************************************************************%

VERSION = '$Revision: 1.6 $';
if ( nargin == 0 ), help(mfilename); return; end;

global DEBUG;   if ( isempty(DEBUG)   ), DEBUG   = 0; end;
global VERBOSE; if ( isempty(VERBOSE) ), VERBOSE = 0; end;

if nargin == 1
    PercentDataUse = 0.5
    RegressionMethod = 2;
elseif nargin == 2
    RegressionMethod = 2;
end
%--------------------------------------------------------------------------%

  % ALGORITHM OUTLINE:
  %
  % (1) average (complex-valued) subset of forward and reverse lines
  % (2) subtract phases of reverse line from forward line FIRST
  % (3) fit a first-order polynomial to the phase difference

  % perform IDFT along readout direction
  projections = mrir_iDFT_freqencode(phascor);

  % crop out extra samples due to oversampling k-space
  projections = mrir_image_crop(projections);
  projections = mrir_image_crop(projections);

  dims = size(projections);

  % ceil is for backward compatibility (if NSeg is 1)
  dims(mrir_DIM_SEG) = ceil( dims(mrir_DIM_SEG)/2 );

  % pick the first several lines (before too much T2* delay)
  % NOTE: this magic number is only a guess
  Nref = min([8, size(projections,2)]);


  % separate out forward and reverse reference lines (assuming polarity
  % reverses every line)

  % NOTE: convention is that segment "0" is a normal line and segment "1"
  % is a reflected line (but it doesn't matter which is which as long as
  % normal and reflected are in different segments)
  %                          1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
  ref_FWD = mean(projections(:,:,:,:,:,:,:,1,:,:,:,:,:,:,:,:), 2);  % line "0" is forward
  ref_REV = mean(projections(:,:,:,:,:,:,:,2,:,:,:,:,:,:,:,:), 2);  % line "1" is reverse

  if ( 0 ),
  ref_o = mean(projections(:,1:2:Nref,:,:,:,:,:,:,:,:,:,:,:,:,:,:), 2);  % line "0" is odd
  ref_e = mean(projections(:,2:2:Nref,:,:,:,:,:,:,:,:,:,:,:,:,:,:), 2);  % line "1" is even
  end;

  % collapse two sets of reference lines each into an Nreadout x Nlines array
  ref_REV = reshape(ref_REV, [dims(1), prod(dims(3:end))]);
  ref_FWD = reshape(ref_FWD, [dims(1), prod(dims(3:end))]);

  % [Nreadout x Nlines]
  phasediff = unwrap( angle(ref_FWD(:,:)./ref_REV(:,:)), [], 1);

  % simple test for 2pi shifts
  % [       1 x Nlines]
  phasediff_center_pos = (mean(phasediff, 1) > +pi) * 2*pi;
  % [Nreadout x Nlines]
  phasediff_offset_pos = repmat(phasediff_center_pos, size(phasediff,1), 1);

  phasediff_center_neg = (mean(phasediff, 1) < -pi) * 2*pi;
  phasediff_offset_neg = repmat(phasediff_center_neg, size(phasediff,1), 1);

  phasediff = phasediff - phasediff_offset_pos;
  phasediff = phasediff + phasediff_offset_neg;
  
  
  %PercentDataUse = 0.6;
  
  % calculate centroid of the data 
  % do it only for first repetition and set of the data!
  intensity = sum(sum(   abs(projections(:,:,:,1,1,1,1,1:2,1,:))   ,8) ,2);
  d = repmat( [1:size(intensity,1)].' , [1, 1, dims(3),1,1,1 1,1,1,dims(10)]) ;
  Centroid = [];
  CenterPoint = ceil( (dims(1)+1)/2 );
  Centroid =  squeeze(  round( sum( intensity(:,1,:,1,1,1,1,1,1,:).*d ,1 ) ./ sum(intensity,1) )   -  CenterPoint  );
  %clear d intensity
 
  PointsUsed = (ceil( size(phasediff,1)*PercentDataUse/2 )*2) -1; % always odd number
  StartPoint = CenterPoint - (PointsUsed-1)/2 ;
  EndPoint = CenterPoint + (PointsUsed-1)/2 ;
  x_refline_centered = [1: PointsUsed].';
  x_refline_centered = x_refline_centered - ceil( (length(x_refline_centered)+1)/2 );
  
  NegEdge = ceil( -floor( dims(1)/2 ) + (PointsUsed+1)*0.5 ); % center point can not go further than this 
  CentroidNegMask = Centroid < NegEdge;
  Centroid(CentroidNegMask) = NegEdge;
  
  PosEdge = floor( floor( (dims(1)-1)/2 ) - (PointsUsed+1)*0.5 ); 
  CentroidPosMask = Centroid > PosEdge;
  Centroid(CentroidPosMask) =  PosEdge; 
  

  linear_fit_coeff = zeros(2,size(phasediff,2));
  counter = 1;
  warning('off','all')
  
  for SlcCount = 1:dims(10) 
      for OthersCount = 1:prod(dims(4:7))
          for CoilCount = 1:dims(3)
              x_refline_current = x_refline_centered + Centroid(CoilCount,SlcCount);
              phasediff_current = phasediff(x_refline_current+CenterPoint,counter);
              linear_fit_coeff(:,counter) = robustfit(x_refline_current,phasediff_current);
              counter = counter + 1;
          end
      end
  end
  
  warning('on','all')
%   
%   for count = 1:32
%       x_ref= x_refline_centered + Centroid(count,8);  
%       figure(50); subplot(5,7,count); plot(phasediff(:,(counter-1)+count));
%       figure(51); subplot(5,7,count); plot(phasediff(x_ref+CenterPoint,(counter-1)+count));
%       figure(52); subplot(5,7,count); plot(intensity(:,1,count,1,1,1,1,1,1,8));
%   end 
%   figure(54); plot(Centroid(:,8));
%     
if(0)
    
    if nargin == 1
        % keep only the center 50% of samples, since noise at ends can bias fit
        % NOTE: this magic number is only a guess
        phasediff = phasediff(ceil(1*end/4):floor(3*end/4),:);
    else
        StartPoint = ceil( ((1-PercentDataUse)/2) *size(phasediff,1));
        EndPoint = StartPoint + floor(PercentDataUse*size(phasediff,1));
        phasediff = phasediff(StartPoint:EndPoint,:); %keep more because of SIR for example
    end
    
    % build indices into readout direction for offset and slope of phasediff
    % (since "phasediff" is truncated, these will not be same size as the data
    % correction terms)
    RegressionMethod = 1;
    x_refline = [ones(size(phasediff,1),1), ([0:size(phasediff,1)-1]-size(phasediff,1)/2)'];
    %x_refline = [([0:size(phasediff,1)-1]-size(phasediff,1)/2)'];
    %% method 1
    
    if RegressionMethod == 1
        
        % compute linear (1st order poly) fit for each reference line
        linear_fit_coeff = pinv(x_refline)*phasediff;
        
        linear_fit_coeff = [linear_fit_coeff(1,:); linear_fit_coeff(2,:)];
        est = x_refline*linear_fit_coeff;
    else
        %% method 2
        warning('off','all')
        linear_fit_coeff = [];
        for CoilCount = 1:size(phasediff,2)
            linear_fit_coeff(:,CoilCount) = robustfit(x_refline(:,2),phasediff(:,CoilCount));
        end
        warning('on','all')
    end
    %% debugging plot
    est = x_refline*linear_fit_coeff;
    
    if(0)
        figure(8); %clf;
        Slice = 8;
        Rep = 3;
        Ncoils = dims(3);
        for CoilCount = 1:Ncoils
            subplot(5,7,CoilCount); plot(phasediff(:,CoilCount + (Rep-1)*Ncoils + (Slice-1)*Ncoils*dims(7)),'k');
            %hold on; plot(est(:,CoilCount + (Rep-1)*Ncoils + (Slice-1)*Ncoils*dims(7)),'b');
            %hold on; plot(est(:,CoilCount + (Rep-1)*Ncoils + (Slice-1)*Ncoils*dims(7)) ,'g');
        end
    end
    
end
%%%%%%%%%
  
  
  
  coeff_avg = mean(linear_fit_coeff, 2);

  if ( VERBOSE || DEBUG ),
    disp(sprintf('<i> [%s]: average parameters:  offset = %2.2f deg; slope = %2.2f deg', ...
                 mfilename, (coeff_avg*180/pi)));
  end;

  if ( DEBUG ),
    disp(sprintf('<d> [%s]: generating debugging figures', mfilename));
    mrir_artifact_ghost_compute__DEBUG(phascor, phasediff, linear_fit_coeff, ref_REV, ref_FWD, dims);
  end;

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_artifact_ghost_compute.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
