function level = otsuthresh_scale(I,scale)
% OTSUTHRESH_SCALE Image threshold using linear or log weighted Otsu's method.
%
% level = otsuthresh_scale(I,scale)
%
% Input
%   I : image. If I is single or double, then it is assumed to be in the
%     range of [0,1].
%   scale : either 'linear' or 'log'. If 'linear', then normal Otsu
%     thresholding will be performed. If 'log', then the histogram centers
%     will be log2 weighted; this returns a threshold similar to that
%     computed by CellProfiler software.
%
% Output
%   level : The computed threshold value

% Mostly taken from Matlab's graythresh and otsuthresh functions.
% James Kapaldo

validateattributes(I,{'uint8','uint16','double','single','int16'},{'nonsparse'}, ...
    mfilename,'I',1);

if ~isempty(I)
    
    % Image class range
    a = getrangefromclass(I);
    if a(2) == 1
        minI = min(I(:));
        maxI = max(I(:));
        I = (I-minI)/(maxI-minI);
    end
    
    % Convert all N-D arrays into a single column.  Convert to uint8 for
    % fastest histogram computation.
    I = im2uint8(I(:));
    num_bins = 256;
    counts = imhist(I,num_bins);
    
    % Make counts a double column vector
    counts = double( counts(:) );
    
    % Variables names are chosen to be similar to the formulas in
    % the Otsu paper.
    p = counts / sum(counts);
    omega = cumsum(p);
    
    if strcmp(scale,'log')
        x = log2((1:num_bins)');
    elseif strcmp(scale,'linear')
        x = (1:num_bins)';
    else
        error('otsuthresh_scale:badInput','scale should be either ''linear'' or ''log''.');
    end
    mu = cumsum(p .* x);
    mu_t = mu(end);
    
    sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));
    
    % Find the location of the maximum value of sigma_b_squared.
    % The maximum may extend over several bins, so average together the
    % locations.  If maxval is NaN, meaning that sigma_b_squared is all NaN,
    % then return 0.
    maxval = max(sigma_b_squared);
    isfinite_maxval = isfinite(maxval);
    if isfinite_maxval
        idx = mean(find(sigma_b_squared == maxval));
        % Normalize the threshold to the range [0, 1].
        level = (idx - 1) / (num_bins - 1);
        % Convert the threshold back into the original image range
        level = a(2) * level;
        
        if a(2) == 1
            % Unscale threshold
            level = level*(maxI-minI) + minI;
        end
    else
        level = 0.0;
    end
else
    level = 0.0;
end
