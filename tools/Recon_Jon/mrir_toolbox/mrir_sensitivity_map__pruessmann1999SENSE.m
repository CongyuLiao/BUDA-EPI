function sens = mrir_sensitivity_map(img, ref, P, Omega, noiselevel)
  

  % don't need to choose range for x,y around a point x0,y0, since the
  % weighting term w(x,y) determines the effective range---as long as
  % support of w(x,y) extents slightly beyond range of (x,y).
  
  
% convert noise level measured from raw data into expected pixel intensities
% in magnitude (i.e., modulus) image
  
  
%  "Phi is determined from the denominator image by pixel-wise modulus
%  discrimination with a threshold on the order of the noise level..."

  
  
% "...and further exclusion of pixels with then sparse neighborhood"