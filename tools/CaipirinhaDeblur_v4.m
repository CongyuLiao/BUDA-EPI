function [Kcorrected] = CaipirinhaDeblur_v4(K, prot, evp , PhaseShiftBtwSimulSlices, SliceSep)
%% ________________________________________________________________________
%% Caipirinha and Deblur for Slice GRAPPA
%% ________________________________________________________________________
%% 
%% Haifeng Wang 
%% Sept. 13 2016
%%
%% now work for all acquisition orientations
%% assume that the slices are already sorted (so can use para in original prot to calc stuff) 
%%
%% SliceSep in mm
%%
%% Sept. 13,2016 by Haifeng Wang
%% Add FOV/8 and FOV/10 shift
%%
%% Mar. 31,2017 by Haifeng Wang
%% Add shifting for CaipirinhaShift_K_v2.m
%% ________________________________________________________________________

% correction of SlicePos value due to rotation of principle gradient axis
if size(K,10) > 1;
    
    % When the value for these coordinate is very small, e.g.
    % sPosition_dSag = 1.343452e-9 then the read_meas_dat would not
    % recognize it and will leave the array empty so fix it here
    if isempty(prot.sSliceArray(1).sPosition_dSag) && isempty(prot.sSliceArray(2).sPosition_dSag)
        for count = 1:length(prot.sSliceArray(1:end))
            prot.sSliceArray(count).sPosition_dSag = 0;
            prot.sSliceArray(count).sNormal_dSag = 0;
        end
    end
    if isempty(prot.sSliceArray(1).sPosition_dCor) && isempty(prot.sSliceArray(2).sPosition_dCor)
        for count = 1:length(prot.sSliceArray(1:end))
            prot.sSliceArray(count).sPosition_dCor = 0;
            prot.sSliceArray(count).sNormal_dCor = 0;
        end
    end
    if isempty(prot.sSliceArray(1).sPosition_dTra) && isempty(prot.sSliceArray(2).sPosition_dTra)
        for count = 1:length(prot.sSliceArray(1:end))
            prot.sSliceArray(count).sPosition_dTra = 0;
            prot.sSliceArray(count).sNormal_dTra = 0;
        end
    end
    
    NormalVec = [prot.sSliceArray(1).sNormal_dSag, prot.sSliceArray(1).sNormal_dCor, prot.sSliceArray(1).sNormal_dTra].';   
    Pos(:,1) = [prot.sSliceArray(1:end).sPosition_dSag].';
    Pos(:,2) = [prot.sSliceArray(1:end).sPosition_dCor].';
    Pos(:,3) = [prot.sSliceArray(1:end).sPosition_dTra].';
    SlicePos = Pos*NormalVec;
else
    keyboard('only single slice data so cant determine correctionFac if Gz is rotated')
end
      
%SlicePos  = SlicePos(end:-1:1); 
PhaseShiftPerMM = (PhaseShiftBtwSimulSlices/SliceSep);

%Kcorrected = single(zeros(size(K)));
Kcorrected = (zeros(size(K)));

if PhaseShiftBtwSimulSlices ~= 0
    for SlcCount = 1: size(K,10)
        if abs(PhaseShiftBtwSimulSlices) == 2*pi/2 % FOV/2 shift
            Kcorrected(:,1:2:end,:,:,:,:,:,:,:,SlcCount) =  K(:,1:2:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,2:2:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:2:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*SlicePos(SlcCount));
            %% NEW! change because of CaipirinhaShift_K_v2
            Kcorrected = Kcorrected*exp(i*(PhaseShiftPerMM*SlicePos(SlcCount))/2); 
         elseif abs(PhaseShiftBtwSimulSlices) == 2*pi/3 % FOV/3 shift
            Kcorrected(:,1:3:end,:,:,:,:,:,:,:,SlcCount) =  K(:,1:3:end,:,:,:,:,:,:,:,SlcCount);    
            Kcorrected(:,2:3:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:3:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,3:3:end,:,:,:,:,:,:,:,SlcCount) =  K(:,3:3:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*2*PhaseShiftPerMM*SlicePos(SlcCount));
        elseif abs(PhaseShiftBtwSimulSlices) == 2*pi/4 % FOV/4 shift
            Kcorrected(:,1:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,1:4:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,2:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:4:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,3:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,3:4:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*2*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,4:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,4:4:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*3*PhaseShiftPerMM*SlicePos(SlcCount));
            %% NEW! change because of CaipirinhaShift_K_v2
            Kcorrected = Kcorrected*exp(i*(PhaseShiftPerMM*SlicePos(SlcCount)*3)/2);
        elseif abs(PhaseShiftBtwSimulSlices) == 2*pi/6 % FOV/6 shift
            Kcorrected(:,1:6:end,:,:,:,:,:,:,:,SlcCount) =  K(:,1:6:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,2:6:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:6:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,3:6:end,:,:,:,:,:,:,:,SlcCount) =  K(:,3:6:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*2*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,4:6:end,:,:,:,:,:,:,:,SlcCount) =  K(:,4:6:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*3*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,5:6:end,:,:,:,:,:,:,:,SlcCount) =  K(:,5:6:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*4*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,6:6:end,:,:,:,:,:,:,:,SlcCount) =  K(:,6:6:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*5*PhaseShiftPerMM*SlicePos(SlcCount));              
        elseif abs(PhaseShiftBtwSimulSlices) == 2*pi/8 % FOV/8 shift
            Kcorrected(:,1:8:end,:,:,:,:,:,:,:,SlcCount) =  K(:,1:8:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,2:8:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:8:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,3:8:end,:,:,:,:,:,:,:,SlcCount) =  K(:,3:8:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*2*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,4:8:end,:,:,:,:,:,:,:,SlcCount) =  K(:,4:8:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*3*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,5:8:end,:,:,:,:,:,:,:,SlcCount) =  K(:,5:8:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*4*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,6:8:end,:,:,:,:,:,:,:,SlcCount) =  K(:,6:8:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*5*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,7:8:end,:,:,:,:,:,:,:,SlcCount) =  K(:,7:8:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*6*PhaseShiftPerMM*SlicePos(SlcCount)); 
            Kcorrected(:,8:8:end,:,:,:,:,:,:,:,SlcCount) =  K(:,8:8:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*7*PhaseShiftPerMM*SlicePos(SlcCount));
            %% NEW! change because of CaipirinhaShift_K_v2
            Kcorrected = Kcorrected*exp(i*(PhaseShiftPerMM*SlicePos(SlcCount)*7)/2);
        elseif abs(PhaseShiftBtwSimulSlices) == 2*pi/10 % FOV/10 shift
            Kcorrected(:,1:10:end,:,:,:,:,:,:,:,SlcCount) =  K(:,1:10:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,2:10:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:10:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,3:10:end,:,:,:,:,:,:,:,SlcCount) =  K(:,3:10:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*2*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,4:10:end,:,:,:,:,:,:,:,SlcCount) =  K(:,4:10:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*3*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,5:10:end,:,:,:,:,:,:,:,SlcCount) =  K(:,5:10:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*4*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,6:10:end,:,:,:,:,:,:,:,SlcCount) =  K(:,6:10:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*5*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,7:10:end,:,:,:,:,:,:,:,SlcCount) =  K(:,7:10:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*6*PhaseShiftPerMM*SlicePos(SlcCount)); 
            Kcorrected(:,8:10:end,:,:,:,:,:,:,:,SlcCount) =  K(:,8:10:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*7*PhaseShiftPerMM*SlicePos(SlcCount));
            Kcorrected(:,9:10:end,:,:,:,:,:,:,:,SlcCount) =  K(:,9:10:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*8*PhaseShiftPerMM*SlicePos(SlcCount)); 
            Kcorrected(:,10:10:end,:,:,:,:,:,:,:,SlcCount) = K(:,10:10:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*9*PhaseShiftPerMM*SlicePos(SlcCount));             
        end
    end
else
    Kcorrected = K;
end
