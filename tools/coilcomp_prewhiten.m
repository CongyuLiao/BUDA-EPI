clear all; close all; clc;

%% Paths ----------------------------------------------------------------------------------------------------
addpath('lib');

%% Parameters -----------------------------------------------------------------------------------------------
data = 'rawdata/meas_MID01001_FID50799_tfl_wshfl_VD_R02.dat';
mask = 'mask/vd_r02';

% Data dimensions.
SX  = 256;                  % Readout.
sx  = 256;                  % Readout after wave-zpad.
SY  = 256;                  % Phase encode 1.
sy  = 256;                  % Phase encode 1 after wave-zpad.
SZ  = 256;                  % Phase encode 1.
sz  = 256;                  % Phase encode 2 after wave-zpad.
cc  = 32;                   % Number of channels pre coil compression.
nc  = 16;                   % Number of channels after coil compression.
tf  = 256;                  % Echo train length.
OS  = 6;                    % Readout over-sampling factor for wave.
OSc = 2;                    % Oversampling factor for GRE reference.
wx  = sx * OS;

% Wave parameters.
cY  = [(-0.006890) + 1i * (-0.591000), ( 0.0000000)];
cZ  = [(-0.598000) + 1i * ( 0.032600), (-1.2000000)];
l   = [9, 1, 9, 1];
dy  = 0.1 * (SY/sy);
dz  = 0.1 * (SZ/sz);
adc = 2.4576;               % ADC duration in ms.
gm  = 1.6;                  % Max gradient amplitude in Gauss/cm.
sm  = 18700;                % Gradient slew rate in Gauss/cm/s.
cyc = 8;                    % Number of wave-cycles.
ofY = -2;                   % Offset in Y. (TODO: Might be -20)
ofZ = 0;                    % Offset in Z.

% Temporal parameters.
esp   = 5.3;                % Echo spacing in ms.
TI    = 1100;               % Inversion time in ms.
ttfa  = TI - (tf/2) * esp;  % Time between 180* and first alpha.
alpha = 8;                  % Flip angle of excitation in degrees.
BW    = 410;                % BW in Hz/Px.

%% Create espirit maps and coil compression matrix ----------------------------------------------------------
if (~isfile('data/MID01001/whiten_cc_mat.cfl'))
  fprintf('Preparing reference and whitening matrix. ---------------------------------------------------\n');
  obj = mapVBVD(data, 'ignoreSeg');

  % Extract reference scan and noise data.
  rref = permute(obj{numel(obj)}.refscan(), [1, 3, 4, 2]);
  nmat = permute(squeeze(obj{numel(obj)}.noise()), [2, 1, 3, 4]);
  nmat = reshape(nmat, size(nmat, 1), prod(size(nmat))/size(nmat, 1));

  % Account for over-sampling in reference scan.
  rref = F_inv(rref, [1, 2, 3]);
  rref = rref(floor((size(rref, 1) * (1 - 1/OS))/2) + [1:size(rref, 1)/OS], :, :, :);
  [rx, ry, rz, rc] = size(rref);

  % SVD to get coil compression matrix.
  rref = reshape(rref, rx * ry * rz, rc).'; 
  [u, s, ~] = svd(rref, 'econ');
  u = u(:, 1:nc);
  s = diag(s);

  % Coil compressing noise.
  nmat = u' * nmat;

  % Estimating whitening matrix.
  covm = (nmat * nmat')/(size(nmat, 2) - 1);
  whmt = inv(chol(covm, 'lower'));

  % Joint coil compression and whitening matrix.
  whtcc = whmt * u';

  % Whiten and coil compress reference data.
  rref = whtcc * rref;
  rref = rref.';
  rref = reshape(rref, rx, ry, rz, nc);
  rref = F_fwd(rref, [1, 2, 3]);

  % Resizing data for BART espirit.
  ref = zeros(sx, sy, sz, nc);
  ref(floor((sx - size(rref, 1))/2) + [1:size(rref, 1)], ...
      floor((sy - size(rref, 2))/2) + [1:size(rref, 2)], ...
      floor((sz - size(rref, 3))/2) + [1:size(rref, 3)], ...
      :) = rref;

  % Saving data.
  writecfl('data/MID01001/reference',     ref);
  writecfl('data/MID01001/whiten_cc_mat', whtcc);
end

%% Temporal basis -------------------------------------------------------------------------------------------
if (~isfile('data/MID01001/mprage_phi.cfl'))
  fprintf('Preparing temporal basis. -------------------------------------------------------------------\n');
  [phi, sig, sv, t1v] = mprage_basis(2, 20:5:5000, ttfa, alpha, esp, tf);
  writecfl('data/MID01001/mprage_phi', phi);
  writecfl('data/MID01001/mprage_sig', sig);
  writecfl('data/MID01001/mprage_sv',  sv);
  writecfl('data/MID01001/mprage_t1v', t1v);
end

%% Temporal basis -------------------------------------------------------------------------------------------
if (~isfile('data/MID01001/mprage_phi_3.cfl'))
  fprintf('Preparing temporal basis. -------------------------------------------------------------------\n');
  [phi, sig, sv, t1v] = mprage_basis(3, 20:5:5000, ttfa, alpha, esp, tf);
  writecfl('data/MID01001/mprage_phi_3', phi);
  writecfl('data/MID01001/mprage_sig_3', sig);
  writecfl('data/MID01001/mprage_sv_3',  sv);
  writecfl('data/MID01001/mprage_t1v_3', t1v);
end

%% Temporal basis -------------------------------------------------------------------------------------------
if (~isfile('data/MID01001/mprage_phi_4.cfl'))
  fprintf('Preparing temporal basis. -------------------------------------------------------------------\n');
  [phi, sig, sv, t1v] = mprage_basis(4, 20:5:5000, ttfa, alpha, esp, tf);
  writecfl('data/MID01001/mprage_phi_4', phi);
  writecfl('data/MID01001/mprage_sig_4', sig);
  writecfl('data/MID01001/mprage_sv_4',  sv);
  writecfl('data/MID01001/mprage_t1v_4', t1v);
end

%% Temporal basis -------------------------------------------------------------------------------------------
if (~isfile('data/MID01001/mprage_phi_5.cfl'))
  fprintf('Preparing temporal basis. -------------------------------------------------------------------\n');
  [phi, sig, sv, t1v] = mprage_basis(5, 20:5:5000, ttfa, alpha, esp, tf);
  writecfl('data/MID01001/mprage_phi_5', phi);
  writecfl('data/MID01001/mprage_sig_5', sig);
  writecfl('data/MID01001/mprage_sv_5',  sv);
  writecfl('data/MID01001/mprage_t1v_5', t1v);
end

%% Load data ------------------------------------------------------------------------------------------------
if (~isfile('data/MID01001/tbl.cfl'))
  fprintf('Preparing ADC table. ------------------------------------------------------------------------\n');
  obj = mapVBVD(data, 'ignoreSeg');
  whtcc = readcfl('data/MID01001/whiten_cc_mat');

  h = read_reorder_header(mask);
  nadc = h(5);

  lin = obj{numel(obj)}.image.Lin - 1 - floor(SY/2) + floor(sy/2);
  par = obj{numel(obj)}.image.Par - 1 - floor(SZ/2) + floor(sz/2);
  eco = mod(0:1:(numel(lin)-1), tf);
  rdr = [lin(:), par(:), eco(:)];
  rdr = rdr(1:nadc, :);

  tbl  = permute(obj{numel(obj)}.image.unsorted(), [2, 1, 3]);
  tbl  = tbl(:, :, 1:nadc);
  rc   = size(tbl, 1);
  nacq = size(tbl, 3);

  tbl = reshape(tbl, rc, wx * nacq);
  tbl = reshape(whtcc * tbl, nc, wx, nacq);
  tbl = permute(tbl, [2, 1, 3]);

  msk = table_to_array(sy, sz, tf, rdr + 1); % Plus one for MATLAB indexing.
  dlt = F_inv(squeeze(sum(msk, 6)), [1, 2]);

  writecfl('data/MID01001/tbl',    tbl);
  writecfl('data/MID01001/rdr',    rdr);
  writecfl('data/MID01001/msk',    msk);
  writecfl('data/MID01001/msksum', squeeze(sum(msk, 6)));
  writecfl('data/MID01001/mskpsf', dlt);

  rdr(:, 3) = 0;
  writecfl('data/MID01001/single_rdr', rdr);
  writecfl('data/MID01001/single_phi',   1);
end

%% Generate PSF ---------------------------------------------------------------------------------------------
if (~isfile('data/MID01001/wave.cfl'))
  fprintf('Preparing wave-psf. -------------------------------------------------------------------------\n');
  c = [real(cY(:)); imag(cY(:)); real(cZ(:)); imag(cZ(:))] * wx;
  y = dy * ((0:sy-1) - floor(sy/2)) - ofY;
  z = dz * ((0:sz-1) - floor(sz/2)) - ofZ;
  [p, py, pz] = wave_coeffs_to_psf(c, l, wx, y, z, 0);
  writecfl('data/MID01001/wave', p);
end

%% Calibrate maps -------------------------------------------------------------------------------------------
if (~isfile('data/MID01001/maps.cfl'))
  fprintf('Please run bart.sh for ESPIRiT maps.\n');
  return;
end

%% Reconstruction -------------------------------------------------------------------------------------------
if (~isfile('res/MID01001/coeffs.cfl'))
  fprintf('Please run bart.sh.\n');
  return;
end

%% ----------------------------------------------------------------------------------------------------------
