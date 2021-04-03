function I_comb = mrir_array_combine_wrap(I,CoilCombMethod,covmtx,receive,NormalizationType)

% Kawin Setsompop
% March 24th 2012
% wrapper around various combinatin methods


if CoilCombMethod == 1
    % RSS
    I_comb= mrir_array_combine(I,0);
elseif CoilCombMethod == 2
    % Sensitivity profile only
    I_comb = ...
        mrir_array_combine_optimalSNR_fast(I, receive , diag(ones(size(I,3),1)) ,NormalizationType);
elseif CoilCombMethod == 3
    % NoiseCov
    I_comb = mrir_array_combine_cov(I, covmtx);
elseif CoilCombMethod == 4
    % Optimal SNR
    I_comb = ...
        mrir_array_combine_optimalSNR_fast(I, receive, covmtx, NormalizationType);
end