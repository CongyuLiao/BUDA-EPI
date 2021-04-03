function img_new= rot_mat(img_old,a)
for nn=1:size(img_old,4)
for ii=1:size(img_old,3)
    if a==90
    img_new(:,:,ii,nn)=rot90(img_old(:,:,ii,nn));
    elseif a==180
    img_new(:,:,ii,nn)=rot90(rot90(img_old(:,:,ii,nn)));
    elseif a==270 || a==-90
    img_new(:,:,ii,nn)=rot90(img_old(:,:,ii,nn),'k');
    elseif a==0
    img_new(:,:,ii,nn)=img_old(:,:,ii,nn);
    else
        disp('please input angle: 0,90,180,270')
    end
end
        
end
end