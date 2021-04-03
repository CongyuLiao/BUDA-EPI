function img_new= rot_mat(img_old,a)
for ii=1:size(img_old,3)
    for nn=1:size(img_old,4)
        for mm=1:size(img_old,5)
    if a==90
    img_new(:,:,ii,nn,mm)=rot90(img_old(:,:,ii,nn,mm));
    elseif a==180
    img_new(:,:,ii,nn,mm)=rot90(rot90(img_old(:,:,ii,nn,mm)));
    elseif a==270 || a==-90
    img_new(:,:,ii,nn,mm)=rot90(img_old(:,:,ii,nn,mm),'k');
    elseif a==0
    img_new(:,:,ii,nn,mm)=img_old(:,:,ii,nn,mm);
    else
        disp('please input angle: 0,90,180,270')
    end
        end
    end
end
end