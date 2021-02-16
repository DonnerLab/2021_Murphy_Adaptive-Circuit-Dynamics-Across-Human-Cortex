function rgb = cie2rgb(x,y,Y)

%X = Y./y.*x;
%Z = Y./y.*(1-x-y);
z = (1-x-y);

rgb = xyz2rgb([x' y' z']).*255;
