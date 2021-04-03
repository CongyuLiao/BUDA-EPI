function [kspace_cor]=coilcomp_prewht(data,whmcc)


s = size(data);
tmp = permute( data , [3 1 2 4 5 6 7 8 9 10]);
tmp = squeeze( tmp);
tmp = reshape ( tmp, s(3),s(1)*s(2)*s(10));
tmp = whmcc * tmp;
ncc = size( whmcc ,1);
tmp = reshape( tmp, ncc,s(1),s(2),s(10));
tmp = permute(tmp, [2 3 1 4]);
kspace_cor(:,:,:,1,1,1,1,1,1,:) = tmp;

end