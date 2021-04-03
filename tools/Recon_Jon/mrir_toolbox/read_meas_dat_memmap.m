
%% 29 October 2009 
% Kawin Setsompop

%% 6 Nov 2009
% shrink the patrefscan_phascor matrix to save room when saving data (line ~ 55)
% fix hard coding of number of phase correct lines assumption (3)
% Kawin Setsompop

% Based on code written by Jennifer McNab and Thomas Witzel on memmap method of
% reading in data quickly 
% read first repetition using read_meas_dat from Jon Polimeni and use memmap
% method to read the rest. 

% six important things: prot,evp,data,data_phascor1d,patrefscan,patrefscan_phascor

%IMPORTANT NOTE: if use Jon's recon chain then remove the deinterleaving
%stuff. 

% the data is arranged in a coil reordered way, except for 7T data which
% will be a bit weird......


function [meas] = read_meas_dat_memmap(filename,ReadFirstRepNeeded_Flag,nRepToRead,BeginRep,ascendingSlice_acq,Data7T)

if nargin == 4
    ascendingSlice_acq = 0;
    Data7T = 0;
end

save_fname = [filename(1:end-4) '_FirstRep.dat'];

if ReadFirstRepNeeded_Flag == 1
    opt.ReturnStruct=1;
    opt.ReadMultipleRepetitions = 0;
    opt.CanonicalReorderCoilChannels = 1; 
    % meas.data is not actually obtained here but the coil reordering flag here will make sure that other stuff (e.g. fft scaling factor etc) will be reorder 
    % The data is reordered later on in this code using header information
    
    disp('read & save First Rep')
    tic
    
    %readin first rep and save prot,evp,(patrefscan,patrefscan_phascor)
    meas_first = read_meas_dat(filename,opt);
    
    if isfield(meas_first, 'patrefscan')
        [meas_first.evp,meas_first.prot,meas_first.patrefscan] = SiemensIPATConventionWorkAround(meas_first.evp,meas_first.prot,meas_first.patrefscan);
        % not sure how data is organize esp. when the number of PE lines is not divisible by the acceleration factor
        % Sol: use Jon's convention by looking at non-zero values in his data
        meas_first.evp.SegOneLines = find(sum(meas_first.data(:,:,1,1,1,1,1,1,1,1),1) ~= 0);
        meas_first.evp.SegTwoLines = find(sum(meas_first.data(:,:,1,1,1,1,1,2,1,1),1) ~= 0);
    else
        meas_first.evp.SegOneLines = 1:2:meas_first.evp.NLinMeas;
        meas_first.evp.SegTwoLines = 2:2:meas_first.evp.NLinMeas;
    end
    meas_first.evp.PhasCorSegSwap = (sum(meas_first.data_phascor1d(:,1,1,1,1,1,1,1,1,1),1) == 0); % for some data set, Jon's code swap the phase correction segment to match data
    meas_first.evp.NPhaseCorLines = size(meas_first.data_phascor1d,2);
    
    if ascendingSlice_acq == 0
        deinterleave = strcmp(meas_first.prot.ucMultiSliceMode, 'MSM_INTERLEAVED');
    else
        meas_first.prot.ucMultiSliceMode = 'MSM_SEQUENTIAL';
        deinterleave = 2;
    end
    if  (isfield(meas_first, 'patrefscan'))
        [datprune, meas_first.patrefscan] = mrir_array_GRAPPA_prune([], meas_first.patrefscan , meas_first.evp); % need to modify Jon's code a bit to use this, o.w. use code below
        %[datprune, meas_first.patrefscan] = mrir_array_GRAPPA_prune(meas.data, meas_first.patrefscan , meas_first.evp);
        clear datprune
        
        %shrink patrefscan_phascor matrix
        IndexSeg1 = find(meas_first.patrefscan_phascor(end/2,:,1,1,1,1,1,1,1,1)); L1 = length(IndexSeg1);
        IndexSeg2 = find(meas_first.patrefscan_phascor(end/2,:,1,1,1,1,1,2,1,1)); L2 = length(IndexSeg2);
        s = size(meas_first.patrefscan_phascor);
        s(2) = L1 + L2;
        patrefscan_phacsorTemp = zeros(s);
        
        [C1 S1] = find(squeeze(meas_first.patrefscan_phascor(end/2,:,1,1,1,1,1,1,1,:)));
        [C2 S2] = find(squeeze(meas_first.patrefscan_phascor(end/2,:,1,1,1,1,1,2,1,:)));
        for SlcCount = 1:s(10)
            patrefscan_phascorTemp(:,IndexSeg1,:,:,:,:,:,1,:,SlcCount) = meas_first.patrefscan_phascor(:,C1(1+(SlcCount-1)*L1:SlcCount*L1),:,:,:,:,:,1,:,SlcCount);
            patrefscan_phascorTemp(:,IndexSeg2,:,:,:,:,:,2,:,SlcCount) = meas_first.patrefscan_phascor(:,C2(1+(SlcCount-1)*L2:SlcCount*L2),:,:,:,:,:,2,:,SlcCount);
        end
        meas_first.patrefscan_phascor = patrefscan_phascorTemp;
        clear patrefscan_phascorTemp C1 C2 S1 S2 L1 L2 IndexSeg1 IndexSeg2 s
        
        if deinterleave == 1
            meas_first.patrefscan = mrir_image_slice_deinterleave(meas_first.patrefscan);
            meas_first.patrefscan_phascor = mrir_image_slice_deinterleave(meas_first.patrefscan_phascor);
        elseif deinterleave == 2 % ascending need to reverse        
            meas_first.patrefscan = meas_first.patrefscan(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
            meas_first.patrefscan_phascor = meas_first.patrefscan_phascor(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
        end
        
    end

    meas_first.data = [];
    meas_first.data_phascor1d = []; 
    
    SaveRawData(meas_first,save_fname);
    
    toc
    
else
    disp('Load First Rep Info')
    tic
    meas_first = ReadRawData(save_fname);
    toc
    if ascendingSlice_acq == 0
        deinterleave = strcmp(meas_first.prot.ucMultiSliceMode, 'MSM_INTERLEAVED');
    else
        meas_first.prot.ucMultiSliceMode = 'MSM_SEQUENTIAL';
        deinterleave = 2;
    end  
end

% Get coil reordering info
[fp] = fopen(filename, 'r', 'l');
data_start = fread(fp, 1, 'uint32');
header = fscanf(fp, '%c', data_start-4);
[coil_index, coil_order] = read_meas_dat__reorder_coil_channels(header);

if(Data7T)
if sum(coil_index - [length(coil_index):-1:1]) == 0
    disp('****************************************')
    disp('****************************************')
    disp('Hack for non-contiguous 7T coil reordering stuff')
    disp('****************************************')
    disp('****************************************')
    coil_index = 1:length(coil_index);
end
end
fclose(fp);
    
%% extract params from first rep

meas = meas_first;
clear meas_first;

sData = [meas.evp.NColMeas, meas.evp.NLinMeas, meas.evp.NChaMeas,...
          1, 1, 1, nRepToRead, 2, 1, meas.evp.NSlcMeas];
sPhaseCor = sData; sPhaseCor(2) = meas.evp.NPhaseCorLines;

nRead = sData(1);
nPE = sData(2);
nCoil = sData(3);
nSlice = sData(10);
nPhaseCor = sPhaseCor(2);

nPE_Raw = meas.evp.RawLin;
nRep = meas.evp.RawRep;
meas.evp.NRepMeas = nRepToRead;

SegOneLines = meas.evp.SegOneLines;
SegTwoLines = meas.evp.SegTwoLines;
PhasCorSegSwap  = meas.evp.PhasCorSegSwap;


%% calculate parameters needed for mmap and readout and preallocate matrices

meas.data = single(zeros(sData));
meas.data_phascor1d = single(zeros(sPhaseCor));

linesperrep = (nPE_Raw+nPhaseCor)*nCoil*nSlice;

%% mmap
disp('Mmap')
tic
m = mmap_mdh_noheader_rep(filename,linesperrep,nRep);
disp(['Time: ' num2str(toc) ' s'])


%% Readout

ICE_RAWDATA_SCALE       = 131072.0;  % 64 ^ 3 / 2
K_ICE_AMPL_SCALE_FACTOR = 80 * 20 * ICE_RAWDATA_SCALE / 65536;
 
for RepCount = 1:nRepToRead
    disp([ 'Reading in Rep:' num2str(BeginRep+RepCount-1) ])
    tic
    k = m.Data(BeginRep+(RepCount-1)).kdata;
    k = reshape(k,[],linesperrep)*K_ICE_AMPL_SCALE_FACTOR;
    k = k(33:2:end,:)+i*k(34:2:end,:); %32 headers
    k = reshape(k, [nRead, nCoil, nPE_Raw+nPhaseCor, nSlice]);
    meas.data(:,SegOneLines,:,:,:,:,RepCount,1,:,:) = permute(k(:,coil_index,nPhaseCor+1:2:end,:),[1 3 2 4]);
    meas.data(:,SegTwoLines,:,:,:,:,RepCount,2,:,:) = permute(k(end:-1:1,coil_index,nPhaseCor+2:2:end,:),[1 3 2 4]); % need to reverse k-line as assume EPI data
    if PhasCorSegSwap == 0
        meas.data_phascor1d(:,1:2:end,:,:,:,:,RepCount,1,:,:) = permute(k(:,coil_index,1:2:nPhaseCor,:),[1 3 2 4]);
        meas.data_phascor1d(:,2:2:end,:,:,:,:,RepCount,2,:,:) = permute(k(end:-1:1,coil_index,2:2:nPhaseCor,:),[1 3 2 4]); % need to reverse k-line as assume EPI data
    else
        meas.data_phascor1d(:,2:2:end,:,:,:,:,RepCount,1,:,:) = permute(k(:,coil_index,2:2:nPhaseCor,:),[1 3 2 4]);
        meas.data_phascor1d(:,1:2:end,:,:,:,:,RepCount,2,:,:) = permute(k(end:-1:1,coil_index,1:2:nPhaseCor,:),[1 3 2 4]); % need to reverse k-line as assume EPI data
    end
    toc
end

disp(' ')
disp(' ')

% if ( deinterleave ),
%     meas.data = mrir_image_slice_deinterleave(meas.data);
%     meas.data_phascor1d = mrir_image_slice_deinterleave(meas.data_phascor1d);
% end

if deinterleave == 1
    meas.data = mrir_image_slice_deinterleave(meas.data);
    meas.data_phascor1d = mrir_image_slice_deinterleave(meas.data_phascor1d);
elseif deinterleave == 2 % ascending need to reverse
    meas.data = meas.data(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
    meas.data_phascor1d = meas.data_phascor1d(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
end

if  (isfield(meas, 'patrefscan'))
    [meas.data] = mrir_array_GRAPPA_prune(meas.data, [] , meas.evp);
end
clear m
    
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

  
  
  
  
  

