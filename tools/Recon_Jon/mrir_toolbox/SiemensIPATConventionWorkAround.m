function [evp,prot,patrefscan] = SiemensIPATConventionWorkAround(evp,prot,patrefscan) 
%% Modify MEAS file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% workaround for unusual Siemens convention #1:
if ( evp.NFirstRefLin == 0 ),
    evp.NFirstRefLin = mrir_ice_dimensions(patrefscan, 'lin') - evp.NRefLin + 1;
end;

% workaround for unusual Siemens convention #2:

% (possibly Siemens fills in with GRAPPA fewer lines than are in FFT, so
% last lines are effectively zero-padded; this could throw off SNR
% calculations, so by overriding this we force "mrir_epi_GRAPPA" to fill
% in same number of k-space lines as there are image lines.)
if ( ~isempty(evp.NAFLin) && (evp.NAFLin == 1) && (prot.ucPhasePartialFourier == 1) && (evp.NLinMeas < evp.NImageLins) ),
    jnotify
    keyboard
    evp.NLinMeas = evp.NImageLins;
end;