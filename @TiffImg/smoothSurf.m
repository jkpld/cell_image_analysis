function Z = smoothSurf(obj,z,type)
% SMOOTHSURF Wraper for smoothSurface() to save typing.

y = obj.yCenters;
x = obj.xCenters;
smoothRadius = obj.Surface_Smoothing_Radius / (obj.blockSize*obj.mmPerPixel);

if nargin > 2
    Z = TiffImg.smoothSurface(x,y,z,smoothRadius,type);
else
    Z = TiffImg.smoothSurface(x,y,z,smoothRadius);
end

%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
