function CHA = ConvexArea(B)
% CONVEXAREA Compute the convex area of an object from its boundary
%
% CHA = ConvexArea(B)
%
% Input
%   B : boundary of object as returned by bwboundaries
%
% Output
%   CHA : convex area of object

% Compute convex hull area
try
    K = convhull(B(:,1),B(:,2));
catch 
    CHA = 0;
    return;
end

hull = B(K,:);

% Approximate the area of the shape in an image with the area of the
% polygon describing the convex hull.
siz = size(hull,1);
CHA = abs(sum( (hull([2:siz, 1],1) - hull(:,1)) .* (hull([2:siz, 1],2) + hull(:,2)))/2); 

end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
