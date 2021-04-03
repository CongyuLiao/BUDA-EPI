function DisplayImage(I,crop)

if nargin == 1
    crop = 1;
end

if crop == 1
    imagesc(flipud(I(end/4:3*end/4,:).'))
else
    imagesc(flipud(I(:,:).'))
end
