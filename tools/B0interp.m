function [B0_map_interp ]= B0interp(B0_map,target_size)
    ori_size=size(B0_map);
    x=linspace(-1,1,ori_size(1));
    y=linspace(-1,1,ori_size(2));
    z=linspace(-1,1,ori_size(3));
    [X,Y,Z] = meshgrid(x,y,z);
    xa = linspace(-1,1,target_size(1));
    ya = linspace(-1,1,target_size(2));
    za = linspace(-1,1,target_size(3));

    [Xa,Ya,Za] = meshgrid(xa,ya,za);
    B0_map_interp = interp3(X,Y,Z,B0_map,Xa,Ya,Za,'cubic');

    clear ii x X xa Xa y Y ya Ya z Z za Za ans

end