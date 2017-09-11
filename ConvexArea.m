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
K = convhull(B(:,1),B(:,2));
hull = B(K,:);

minCorner = min(B);
maxCorner = max(B);
Bbox = [minCorner, (maxCorner-minCorner) + 1];

y = hull(:,2) - Bbox(2) + 1;
x = hull(:,1) - Bbox(1) + 1;
M = Bbox(4);
N = Bbox(3);

convexImage = poly2mask(x,y,M,N);
CHA = sum(convexImage(:));

end