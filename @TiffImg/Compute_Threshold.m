function th = Compute_Threshold(tiffImg, removeBackgroundFirst, removeForegroundFirst)
% COMPUTE_THRESHOLD Compute the image threshold in each block of the image.
%
% threshold = Compute_Threshold(tiffImg, removeBackgroundFirst, removeForegroundFirst)
%
% Input
%   removeBackgroundFirst : logical flag. If true, and if the TiffImg has a
%     Background, then the background will be removed (subtracted away)
%     before computing the image threshold. (Default, false)
%   removeForegroundFirst : logical flag. If true, and if the TiffImg has a
%     Foreground, then the foreground will be removed (divided out) before
%     computing the image threshold. (Default, false)
%
% Output
%   threshold : Matrix giving the smoothed threshold computed accross the
%     image. Even if no output is requested the threshold is still stored
%     to the tiffImg tiffImgect.
%
% Note : If tiffImg.Surface_Smoothing_Radius is non-NaN, then the threshold
% will be smoothed with a Lowess smoothing surface.

% James Kapaldo

if nargin < 2
    removeBackgroundFirst = false;
end
if nargin < 3
    removeForegroundFirst = false;
end

try
    threshold = zeros(tiffImg.numBlcks,tiffImg.workingClass);
    removeBackgroundFirst = removeBackgroundFirst && ~isempty(tiffImg.BG_smooth);
    removeForegroundFirst = removeForegroundFirst && ~isempty(tiffImg.FG_smooth);
    if removeBackgroundFirst
        BG_fun = generateFunction(tiffImg, tiffImg.BG_smooth, tiffImg.BG_Xstripe, true);
    end
    if removeForegroundFirst
        FG_fun = generateFunction(tiffImg, tiffImg.FG_smooth, tiffImg.FG_Xstripe, true);
    end
    
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose,'name','Computing threshold,');
    progress.start();
    
    % Iterage over x blocks
    for blck_x = 1:tiffImg.numBlcks(2)
        
        tiffImg.open() % Open Image
        
        % Iterate over y blocks
        for blck_y = 1:tiffImg.numBlcks(1)
            
            % Read in image block
            [I,x,y] = getBlock(tiffImg, blck_x, blck_y);

            % Filter image
            if tiffImg.Use_GPU
                I = gpuArray(I);
                Is = imfilter(I,tiffImg.Image_Smooth_Kernel,'symmetric'); clear I
                I = gather(Is); clear Is;
            else
                I = imfilter(I,tiffImg.Image_Smooth_Kernel,'symmetric');
            end
            
            % Remove background
            if removeBackgroundFirst
                I = I - BG_fun(x,y);
            end
            if removeForegroundFirst
                I = I ./ FG_fun(x,y);
            end

            % Compute threshold
%             if any(I(:)>1)
%                 fprintf('larger than 1!\n')
%             end
            threshold(blck_y,blck_x) = tiffImg.otsuthresh_scale(I,'log');
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        progress.iteration_end(); % Update progress counter
    end % x block
    
    % Increase the threshold slightly. The log scale is a bit low.
    threshold = 1.2 * threshold;
    
    % Smooth the threshold surface
    threshold = smoothSurf(tiffImg, threshold);
    
    % Save output for later reference
    tiffImg.threshold = struct('x',tiffImg.xCenters,'y',tiffImg.yCenters,'Z',threshold,'removeBackgroundFirst',removeBackgroundFirst);
    
    % Create function for internally evaluating the threshold in other
    % functions.
    tiffImg.threshold_fun = generateFunction(tiffImg, threshold);
    
    % Modify a flags letting other functions know if they need to apply the
    % threshold after the background correction and foreground correction.
    tiffImg.Threshold_After_Background = removeBackgroundFirst;
    tiffImg.Threshold_After_Foreground = removeForegroundFirst;
    
    % Output the threshold matrix if requested.
    if nargout > 0
        th = threshold;
    end
    
catch ME
    tiffImg.close();
    if tiffImg.Use_GPU
        gpuDevice([]);
    end
    fprintf('entered catch statement\n')
    rethrow(ME)
end


tiffImg.close();

end % function
