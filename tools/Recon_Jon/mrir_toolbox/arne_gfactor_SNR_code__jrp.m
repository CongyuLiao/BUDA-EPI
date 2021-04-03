%--------------------------------------------------------------------------

% Subfunction mang - line-by-line SENSE SNR evaluation

%--------------------------------------------------------------------------





function[SNR_LOSS_line,SNR_SENSE_line,SNR_FULL_line,SNR_ART_line] = mang(b,ref,R,noisecor,corg1,corg2);

%

%   Line by line calculation of SENSE SNR and Full Array SNR

%

%   Is called by: gfcalc.m

%

% Usage[SNR_LOSS_line,SNR_SENSE_line,SNR_FULL_line,SNR_ART_line] = mang(b,ref,R,noisecor,corg1,corg2)

%

% Input variables:

% b             Single line of signal data from all channels

% ref           Single line of sensitivity data (reference) from all channels

% R          Total reduction in both directions: R = red_rd*red_ph

% noisecor      Noise correlation matrix

% corg1:        Conditioning factor for m-SENSE

% corg2:        Conditioning factor for final SNR Loss Map

%

% Output variable(s):

% SNR_LOSS_line      Single line of Relative SNR loss = 1/g-factor

% SNR_SENSE_line     Single line of SENSE SNR

% SNR_FULL_line      Single line of Full Array SNR

% SNR_ART_line       Single line of Artifact Power

%

% Last change:

%   When            What            Who

% 25.08.05      add comments    Arne Reykowski





[nl,nCh] = size(b);           % nCh = number of channels

% nl = Full matrix size in phase encoding

% direction



red_FOV = round(nl/R);   %  Reduced matrix size in phase encoding direction



b   = reshape(b,   red_FOV, R, nCh);  % Separate images into nCh partitions

ref = reshape(ref, red_FOV, R, nCh);  % Separate images into nCh partitions



p = permute(b,[3 2 1]);

p_ref = permute(ref, [3 2 1]);



reg = eye(R)*corg1;   % Matrix for regularization

full_snr_opt = zeros(red_FOV,R);

full_snr_ref = zeros(red_FOV,R);

g_fact_ref = zeros(red_FOV,R);

g_fact = zeros(red_FOV,R);

with_artefact = zeros(red_FOV,R);

sense_snr_ref = zeros(red_FOV,R);



ncor_inv = noisecor^-1;



for i3 = 1:red_FOV

    p1 = p(:,:,i3);                              % Image data

    p1_ref = p_ref(:,:,i3);                      % Reference data

    %-----------------------------------------

    A = p1_ref'*ncor_inv;

    B = A*p1_ref;

    weight0 = A;                                % weighting factor with reference image

    w_norm0 = diag(diag(sqrt(B)).^-1);             % Normalize for weight

    weight1 = w_norm0*weight0;                       %

    ss1 = weight1*p1;                            % Combined signal using reference

    full_snr_ref(i3,:) = abs(diag(ss1).');            % Full SNR with reference weighting (from ss1)

    %---------------------------------------

    clear w_norm0
    clear weight1
    
    Up = ((B+reg)^-1)*(A);

    Ak = diag(diag(sqrt(Up*noisecor*Up')).^-1);

    Uk = Ak*Up;

    rho = Uk*p1;

    sense_snr_ref(i3,:) = abs(diag(rho).');            % accelerated SNR (from ss2)

    with_artefact(i3,:) = abs(sum(rho,2))';

    %-----------------------------------------



end



SNR_FULL_line = abs(reshape(full_snr_ref,nl,1));

SNR_SENSE_line = abs(reshape(sense_snr_ref,nl,1));

reg2 = ones(size(SNR_FULL_line)).*((corg2)^2);

SNR_LOSS_line = sqrt((SNR_SENSE_line.^2+reg2)./(SNR_FULL_line.^2+reg2));

SNR_ART_line = abs(reshape(with_artefact,nl,1)-SNR_SENSE_line);
