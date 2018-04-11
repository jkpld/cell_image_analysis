% INTERP2MEX_WEXPAND Fast 2-D bilinear interpolation with grid expansion
%
% Zi = interp2mex_wExpand(Z, xi, yi)
%
% Interpolates to find Zi (n x m, matrix), the values for the underlying
% 2-D function Z at the grid of points formed by vectors xi (1 x m) and yi
% (n x 1). Out of range values are interpolated to nearest neighbor. Only
% input matrices of class double are supported. Other classes will give
% unpredictable results.
%
% This function assumes that the spacing between all points is 1.
%
% The code is built starting from the function posted on the fileexchange
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/68708
% by Jøger Hansegård

% James Kapaldo





% If no compiler is available to compile interp2mex_wExpand.c, then we
% would need to include checks everywhere in the code to determine which
% interpolation function to use. Instead of doing that, here, we simply
% overload the interp2mex function. If the .mex function is available, then
% it will be called; however, if it is not available, then the below code
% will run.
function zi = interp2mex_wExpand(Z, xi, yi)
% 2D linear interpolation with linear extrapolation.
% Note that the interp2mex_wExpand.mex function uses nearest neighbor
% extrapolation instead of linear.
sz = size(Z);
dx = {1:sz(1), 1:sz(2)};

Z = griddedInterpolant(dx,Z);
zi = Z({yi,xi});

end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
