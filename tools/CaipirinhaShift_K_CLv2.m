function [Kcorrected] = CaipirinhaShift_K_CLv2(K, CurrentSliceGroup,PhaseShiftBase)


PhaseShift = rem(PhaseShiftBase*(CurrentSliceGroup-1), 2*pi);
% NOTE: shift convention follow from deblurring where deblur is done in realtion to the first slice group (ie. first slc group is centered).
%       This is due to the fact that the slicePos param of SMS data is relative to the first slice group (this convention is created in memmap)


Kcorrected = (zeros(size(K)));

if PhaseShift ~= 0
    %     if abs(PhaseShiftBase) == pi % FOV/2 shift
    %         Kcorrected(:,1:2:end,:,:,:,:,:,:,:,:) =  K(:,1:2:end,:,:,:,:,:,:,:,:);
    %         Kcorrected(:,2:2:end,:,:,:,:,:,:,:,:) =  K(:,2:2:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
    %     elseif abs(PhaseShiftBase) == 2*pi/3 % FOV/3 shift
    %         Kcorrected(:,1:3:end,:,:,:,:,:,:,:,:) =  K(:,1:3:end,:,:,:,:,:,:,:,:);
    %         Kcorrected(:,2:3:end,:,:,:,:,:,:,:,:) =  K(:,2:3:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
    %         Kcorrected(:,3:3:end,:,:,:,:,:,:,:,:) =  K(:,3:3:end,:,:,:,:,:,:,:,:)*exp(i*2*PhaseShift);
    %     elseif abs(PhaseShiftBase) == pi/2 % FOV/4 shift
    %         Kcorrected(:,1:4:end,:,:,:,:,:,:,:,:) =  K(:,1:4:end,:,:,:,:,:,:,:,:);
    %         Kcorrected(:,2:4:end,:,:,:,:,:,:,:,:) =  K(:,2:4:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
    %         Kcorrected(:,3:4:end,:,:,:,:,:,:,:,:) =  K(:,3:4:end,:,:,:,:,:,:,:,:)*exp(i*2*PhaseShift);
    %         Kcorrected(:,4:4:end,:,:,:,:,:,:,:,:) =  K(:,4:4:end,:,:,:,:,:,:,:,:)*exp(i*3*PhaseShift);
    %     elseif abs(PhaseShiftBase) == pi/3 % FOV/6 shift
    %         Kcorrected(:,1:6:end,:,:,:,:,:,:,:,:) =  K(:,1:6:end,:,:,:,:,:,:,:,:);
    %         Kcorrected(:,2:6:end,:,:,:,:,:,:,:,:) =  K(:,2:6:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
    %         Kcorrected(:,3:6:end,:,:,:,:,:,:,:,:) =  K(:,3:6:end,:,:,:,:,:,:,:,:)*exp(i*2*PhaseShift);
    %         Kcorrected(:,4:6:end,:,:,:,:,:,:,:,:) =  K(:,4:6:end,:,:,:,:,:,:,:,:)*exp(i*3*PhaseShift);
    %         Kcorrected(:,5:6:end,:,:,:,:,:,:,:,:) =  K(:,5:6:end,:,:,:,:,:,:,:,:)*exp(i*4*PhaseShift);
    %         Kcorrected(:,6:6:end,:,:,:,:,:,:,:,:) =  K(:,6:6:end,:,:,:,:,:,:,:,:)*exp(i*5*PhaseShift);
    %     elseif abs(PhaseShiftBase) == pi/4 % FOV/8 shift
    %         Kcorrected(:,1:8:end,:,:,:,:,:,:,:,:) =  K(:,1:8:end,:,:,:,:,:,:,:,:);
    %         Kcorrected(:,2:8:end,:,:,:,:,:,:,:,:) =  K(:,2:8:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
    %         Kcorrected(:,3:8:end,:,:,:,:,:,:,:,:) =  K(:,3:8:end,:,:,:,:,:,:,:,:)*exp(i*2*PhaseShift);
    %         Kcorrected(:,4:8:end,:,:,:,:,:,:,:,:) =  K(:,4:8:end,:,:,:,:,:,:,:,:)*exp(i*3*PhaseShift);
    %         Kcorrected(:,5:8:end,:,:,:,:,:,:,:,:) =  K(:,5:8:end,:,:,:,:,:,:,:,:)*exp(i*4*PhaseShift);
    %         Kcorrected(:,6:8:end,:,:,:,:,:,:,:,:) =  K(:,6:8:end,:,:,:,:,:,:,:,:)*exp(i*5*PhaseShift);
    %         Kcorrected(:,7:8:end,:,:,:,:,:,:,:,:) =  K(:,7:8:end,:,:,:,:,:,:,:,:)*exp(i*6*PhaseShift);
    %         Kcorrected(:,8:8:end,:,:,:,:,:,:,:,:) =
    %         K(:,8:8:end,:,:,:,:,:,:,:,:)*exp(i*7*PhaseShift);
    %     elseif abs(PhaseShiftBase) == pi/6 % FOV/12 shift
    %         Kcorrected(:,1:12:end,:,:,:,:,:,:,:,:) =  K(:,1:12:end,:,:,:,:,:,:,:,:);
    %         Kcorrected(:,2:12:end,:,:,:,:,:,:,:,:) =  K(:,2:12:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
    %         Kcorrected(:,3:12:end,:,:,:,:,:,:,:,:) =  K(:,3:12:end,:,:,:,:,:,:,:,:)*exp(i*2*PhaseShift);
    %         Kcorrected(:,4:12:end,:,:,:,:,:,:,:,:) =  K(:,4:12:end,:,:,:,:,:,:,:,:)*exp(i*3*PhaseShift);
    %         Kcorrected(:,5:12:end,:,:,:,:,:,:,:,:) =  K(:,5:12:end,:,:,:,:,:,:,:,:)*exp(i*4*PhaseShift);
    %         Kcorrected(:,6:12:end,:,:,:,:,:,:,:,:) =
    %         K(:,6:12:end,:,:,:,:,:,:,:,:)*exp(i*5*PhaseShift);
    %         Kcorrected(:,7:12:end,:,:,:,:,:,:,:,:) =  K(:,7:12:end,:,:,:,:,:,:,:,:)*exp(i*6*PhaseShift);
    %         Kcorrected(:,8:12:end,:,:,:,:,:,:,:,:) =  K(:,8:12:end,:,:,:,:,:,:,:,:)*exp(i*7*PhaseShift);
    %         Kcorrected(:,9:12:end,:,:,:,:,:,:,:,:) =
    %         K(:,9:12:end,:,:,:,:,:,:,:,:)*exp(i*8*PhaseShift);
    %         Kcorrected(:,10:12:end,:,:,:,:,:,:,:,:) =  K(:,10:12:end,:,:,:,:,:,:,:,:)*exp(i*9*PhaseShift);
    %         Kcorrected(:,11:12:end,:,:,:,:,:,:,:,:) =  K(:,11:12:end,:,:,:,:,:,:,:,:)*exp(i*10*PhaseShift);
    %         Kcorrected(:,12:12:end,:,:,:,:,:,:,:,:) =  K(:,12:12:end,:,:,:,:,:,:,:,:)*exp(i*11*PhaseShift);
    %
    %
    Steps = 2*pi/abs(PhaseShiftBase); 
         for count = 1:Steps
             Kcorrected(:,count:Steps:end,:,:,:,:,:,:,:,:) =  K(:,count:Steps:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift*(count-1));
         end
    %         elseif abs(PhaseShift) == 4*pi/3 % FOV/3 shift
    %             Kcorrected(:,1:3:end,:,:,:,:,:,:,:,:) =  K(:,1:3:end,:,:,:,:,:,:,:,:);
    %             Kcorrected(:,2:3:end,:,:,:,:,:,:,:,:) =  K(:,2:3:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
    %             Kcorrected(:,3:3:end,:,:,:,:,:,:,:,:) =  K(:,3:3:end,:,:,:,:,:,:,:,:)*exp(i*2*PhaseShift);
    %   else
    %       keyboard
    %   end
else
    Kcorrected = K;
end

