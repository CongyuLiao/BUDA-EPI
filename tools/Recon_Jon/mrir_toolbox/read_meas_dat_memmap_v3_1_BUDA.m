
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

%% 6 July 2011
% automatically detect if this file has been read before (i.e. no need for
% the ReadFistRepNeeded_Flag)

%% 19 Jan 2012
% change the way memmap work by using pointer associate with data_begin rather than end of file
% this will allow for reconstuct of incomplete files 
%
% make it work for SMS data 
% make it work for both VB17 and VD11
% NOTE: ReorderChannels option is not included here so will not work with old data that have already got .mat saved of first rep!!!

% using Jon's new read_meas_dat that arrange things in cell array (much faster read)

function [meas] = read_meas_dat_memmap_v3_1_BUDA(filename,nRepToRead,BeginRep,SMSdata,ascendingSlice_acq,data_isvd,AccZ,evp)

% if nargin == 4
%     ascendingSlice_acq = 0; 
% end
% if nargin<7
%     AccZ=3;
% end
save_fname = [filename(1:end-4) '_FirstRep.dat'];

if exist([save_fname(1:end-4) '_Raw.mat'],'file') ~= 2
    opt.ReturnStruct=1;
    opt.ReadMultipleRepetitions = 0;
    opt.SqueezeChannels = 1;
    opt.ReturnCellArray = 0;
    opt.DataISVD = data_isvd;
    
    % grab meas.data as a cell which make things a lot faster but some cells can be empty due to PAT and SMS
    if SMSdata == 1
        opt.SMSRefscan = 1;
    else 
        opt.SMSRefscan = 0;
    end
    
    disp('read & save First Rep')
    tic
     
    %readin first rep and save prot,evp,(patrefscan,patrefscan_phascor)
    meas_first = read_meas_dat(filename,opt);
    
    disp('Hack for VD for now!!!*********************************************')
%     CL: overwrite the evp file
    meas_first.evp = evp;
    meas_first.evp.NChaMeas = size(meas_first.data,3);
    disp('*********************************************')
    
    
    if isfield(meas_first, 'patrefscan')
        [meas_first.evp,meas_first.prot,meas_first.patrefscan] = SiemensIPATConventionWorkAround(meas_first.evp,meas_first.prot,meas_first.patrefscan);
        % not sure how data is organize esp. when the number of PE lines is not divisible by the acceleration factor
        % Sol: use Jon's convention by looking at non-zero values in his data
        if( opt.ReturnCellArray == 0)
            % use this if meas.data is a matrix (In Jon's new fast read_meas_dat meas.data is a cell)
            meas_first.evp.SegOneLines = find(sum(meas_first.data(:,:,1,1,1,1,1,1,1,1),1) ~= 0);
            meas_first.evp.SegTwoLines = find(sum(meas_first.data(:,:,1,1,1,1,1,2,1,1),1) ~= 0);
        else
            % use this if meas.data is a cell
            if SMSdata == 1 % find non empty slice and grab the first SlcGroup slc position
                [SlcMask] = SlcMaskGenerator(meas_first);
                index = find(SlcMask == 1);
                FirstIndex = index(1);
                
                Nslices = size(meas_first.smsrefscan,10);
                %NslicesEX =  meas_first.prot.sWiPMemBlock_adFree(1);
                %SlcsPerGroup = Nslices/NslicesEX;
                SlcsPerGroup = sum(SlcMask);
                NslicesEX = Nslices/SlcsPerGroup;
                meas_first.prot.sSliceArray = meas_first.prot.sSliceArray(1:SlcsPerGroup); % CONVENTION: slices position of the top slice group!!!
            else
                FirstIndex = 1;
            end
            b = meas_first.data(1,:,1,1,1,1,1,1,1,FirstIndex);
            b = squeeze(b);
            SegOneLinesMask = length(b);
            for count = 1:length(b)
                SegOneLinesMask(count) = ~isempty(b{count});
            end
            meas_first.evp.SegOneLines = find(SegOneLinesMask == 1);
            
            c = meas_first.data(1,:,1,1,1,1,1,2,1,FirstIndex);
            c = squeeze(c);
            SegTwoLinesMask = length(c);
            for count = 1:length(c)
                SegTwoLinesMask(count) = ~isempty(c{count});
            end
            meas_first.evp.SegTwoLines = find(SegTwoLinesMask == 1);
        end
    else
        meas_first.evp.SegOneLines = 1:2:meas_first.evp.NLinMeas;
        meas_first.evp.SegTwoLines = 2:2:meas_first.evp.NLinMeas;    opt.ReturnStruct = 1;
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
    %    [datprune, meas_first.patrefscan] = mrir_array_GRAPPA_prune(meas_first.data, meas_first.patrefscan , meas_first.evp);
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
    end
    
    if deinterleave == 1
        if SMSdata == 1
            meas_first.smsrefscan = mrir_image_slice_deinterleave(meas_first.smsrefscan);
            meas_first.smsrefscan_phascor = mrir_image_slice_deinterleave(meas_first.smsrefscan_phascor);
        end
        if  (isfield(meas_first, 'patrefscan'))
            meas_first.patrefscan = mrir_image_slice_deinterleave(meas_first.patrefscan);
            meas_first.patrefscan_phascor = mrir_image_slice_deinterleave(meas_first.patrefscan_phascor);
        end
    elseif deinterleave == 2 % ascending need to reverse
        if SMSdata == 1
            meas_first.smsrefscan = meas_first.smsrefscan(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
            meas_first.smsrefscan_phascor = meas_first.smsrefscan_phascor(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
        end
        if  (isfield(meas_first, 'patrefscan'))
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

coil_index = 1:meas_first.evp.NChaMeas;

%% extract params from first rep

meas = meas_first;
clear meas_first;
if(0)
    if SMSdata == 1
        if meas.isvd == 1
            if ~isempty(meas.prot.sWiPMemBlock_adFree)
                meas.prot.sWiPMemBlock_adFree(1) = meas.prot.sWipMemBlock_adFree(5) % SMS factor
                meas.prot.sWiPMemBlock_adFree(2) = meas.prot.sWipMemBlock_adFree(3)*meas.prot.lAccelFactPE % FOV shift factor
            else
                disp('cant find SMS and shift paras in .dat')
                disp('manual entering here.....')
                keyboard
                meas.prot.sWiPMemBlock_adFree(1) = 3;
                meas.prot.sWiPMemBlock_adFree(2) = 3*2;
            end
            meas.prot.sWiPMemBlock_adFree(3) = meas.prot.dThickness/meas.prot.sWiPMemBlock_adFree(1); % slice seperation in mm
            meas.prot.sWiPMemBlock_adFree
        else
            meas.prot.sWiPMemBlock_adFree(1) = meas.prot.sWiPMemBlock_adFree(2); % SMS factor
            meas.prot.sWiPMemBlock_adFree(2) = meas.prot.sWiPMemBlock_adFree(3)*meas.prot.lAccelFactPE; % FOV shift factor
            meas.prot.sWiPMemBlock_adFree(3) = meas.prot.dThickness/meas.prot.sWiPMemBlock_adFree(1); % slice seperation in mm
            meas.prot.sWiPMemBlock_adFree
        end
        NslicesEX =  meas.prot.sWiPMemBlock_adFree(1);
    else
        NslicesEX = 1;
    end
end

NslicesEX = AccZ;

sData = [meas.evp.NColMeas, meas.evp.NLinMeas, meas.evp.NChaMeas,...
          1, 1, 1, nRepToRead, 2, 1, meas.evp.NSlcMeas/NslicesEX ];
sPhaseCor = sData; sPhaseCor(2) = meas.evp.NPhaseCorLines;

nRead = sData(1);
nPE = sData(2);
nCoil = sData(3);
nSlice = sData(10);
nPhaseCor = sPhaseCor(2);

nPE_Raw = meas.evp.RawLin*2-1;
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
%m = mmap_mdh_noheader_rep_v2(filename,linesperrep,nRep,nRead,nCoil,meas.data_begin,meas.isvd);
m = mmap_mdh_noheader_rep_v2(filename,linesperrep,nRep,nRead,nCoil,meas.data_begin,data_isvd);
% note: meas.isvd is added into read_meas_dat by Kawin to have flag to see if vd or vb
%m = mmap_mdh_noheader_rep(filename,linesperrep,nRep);
disp(['Time: ' num2str(toc) ' s'])


%% Readout

ICE_RAWDATA_SCALE       = 131072.0;  % 64 ^ 3 / 2
K_ICE_AMPL_SCALE_FACTOR = 80 * 20 * ICE_RAWDATA_SCALE / 65536;
%if meas.isvd == 1
if data_isvd == 1    
    mdh_length_float32 = 48; % per PE line 
    mdh_ch_length_float32 = 8; % per chanel in each PE line
else
    mdh_length_float32 = 0;
    mdh_ch_length_float32 = 32;
end


for RepCount = 1:nRepToRead
    disp([ 'Reading in Rep:' num2str(BeginRep+RepCount-1) ])
    tic
    k = m.Data(BeginRep+(RepCount-1)).kdata;
%     k = reshape(k,[],linesperrep)*K_ICE_AMPL_SCALE_FACTOR;
%     k = k(33:2:end,:)+i*k(34:2:end,:); %32 headers
%     k = reshape(k, [nRead, nCoil, nPE_Raw+nPhaseCor, nSlice]);
    k = reshape(k,[],linesperrep/nCoil)*K_ICE_AMPL_SCALE_FACTOR;   
    k = k(1+mdh_length_float32:end,:);
    k = reshape(k,[],linesperrep);    
    k = k(1+mdh_ch_length_float32:2:end,:) + i*k(2+mdh_ch_length_float32:2:end,:); %32 headers
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

% if  (isfield(meas, 'patrefscan'))
%     [meas.data] = mrir_array_GRAPPA_prune(meas.data, [] , meas.evp);
% end
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

function [SlcMask] = SlcMaskGenerator(meas_first)

% figure out SMS factor and slc location

if  isfield(meas_first, 'patrefscan')
    R_inplane = meas_first.prot.lAccelFactPE;
else
    R_inplane = 1;
end

PElineSearchCount = 1; 
while PElineSearchCount <= R_inplane+1
    a = meas_first.data(1,PElineSearchCount,1,1,1,1,1,1,1,:);
    a = squeeze(a);
    SlcMask = length(a);
    for count = 1:length(a)
        SlcMask(count) = ~isempty(a{count});
    end
    if sum(SlcMask) ~= 0
        break; % found data!!
    else
        if (PElineSearchCount == R_inplane+1)
            disp('error!! cant find the data')
            keyboard
        end
        PElineSearchCount = PElineSearchCount+ 1; % keep looking!!
    end
end

%SMS_factor = size(meas_first.smsrefscan,10)/sum(SlcMask);
    
    
    
    
    
     
  
  
  
  

