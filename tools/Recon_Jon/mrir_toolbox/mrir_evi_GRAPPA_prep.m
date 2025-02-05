function [dat_corr_grid, varargout] = mrir_evi_GRAPPA_prep(varargin)
%MRIR_EVI_GRAPPA_PREP  perform phase correction and regrid if ramp sampling
%
% [dat_prep, acs_prep] = mrir_evi_GRAPPA_prep(dat, dat_phascor, acs, acs_phascor, prot, evp);

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/sep/06
% $Id: mrir_evi_GRAPPA_prep.m,v 1.3 2008/04/07 00:14:40 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( isstruct(varargin{1}) ),
    meas = varargin{1};

    dat		 = meas.data;
    dat_phascor	 = meas.data_phascor1d;
    acs 	 = meas.patrefscan;
    acs_phascor	 = meas.patrefscan_phascor;
    prot	 = meas.prot;
    evp		 = meas.evp;
    
  else,

    dat		 = varargin{1}; 
    dat_phascor	 = varargin{2};
    acs 	 = varargin{3};
    acs_phascor	 = varargin{4};
    prot	 = varargin{5};
    evp		 = varargin{6};
    
  end;
  
  
  %==--------------------------------------------------------------------==%

  [dat_raw, acs_raw] = mrir_array_GRAPPA_prune(dat, acs, evp);

  %%%  dat_raw = mrir_evi_distort_B0_phaseref(dat_raw, dat_phascor_par);

  dat_corr_grid = mrir_evi_GRAPPA_prep__phascor_regrid(dat_raw, dat_phascor, prot);

  if ( (nargout >= 2) && (~isempty(acs)) ),
    acs_corr_grid = mrir_evi_GRAPPA_prep__phascor_regrid(acs_raw, acs_phascor, prot);
    varargout{1} = acs_corr_grid;
  end;

  return;



%**************************************************************************%
function raw_grid = mrir_evi_GRAPPA_prep__phascor_regrid(raw, phascor, prot)


  % correct EPI ghosting:
  phascor_coeff = mrir_artifact_ghost_compute(phascor);

  hyb_roft = mrir_iDFT_freqencode(raw);
  hyb_corr = mrir_artifact_ghost_correct(hyb_roft, phascor_coeff);

  hyb_coll = mrir_multishot_segment_collapse(hyb_corr);


  %correct ramp sampling by regridding:

  %skip the regridding step if neither trapezoid or sinusoid
  if (prot.alRegridMode == 1 ),
      hyb_grid = hyb_coll;
  else,
      hyb_grid = mrir_regrid_trapezoid(hyb_coll, prot);
  end;

  raw_grid = mrir_fDFT_freqencode(hyb_grid);



  %%% DEBUG:
  if ( 0 ),
    hyb_corr_peft = mrir_iDFT_phasencode(hyb_coll, 'lin');
    img_corr_paft = mrir_iDFT_phasencode(hyb_corr_peft, 'par');
    img_corr_rss  = mrir_array_combine_rss(img_corr_paft);
    
    
    hyb_grid_peft = mrir_iDFT_phasencode(hyb_grid, 'lin');
    img_grid_paft = mrir_iDFT_phasencode(hyb_grid_peft, 'par');
    img_grid_rss  = mrir_array_combine_rss(img_grid_paft);
  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_evi_GRAPPA_prep.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: