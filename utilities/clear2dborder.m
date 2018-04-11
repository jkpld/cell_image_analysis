function bw_clear = clear2dborder(bw)
% CLEAR2DBORDER Fast imclearborder for 2d mask (logical) images. 
%
% bw_clear = clear2dborder(bw)
%
% Input
%   bw : input logical image (mask). 
%
% Output 
%   bw_clear : input logical image with all components touching the border
%     removed.
%
% The input mask can be a gpuArray. If it is a gpuArray, then bw_clear will
% also be a gpuArray.

% James Kapaldo
sz = size(bw);

marker = bw;
marker(2:sz(1)-1,2:sz(2)-1) = false;

on_border = imreconstruct(marker, bw);
bw_clear = bw & ~on_border;
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
