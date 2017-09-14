function th = Compute_Threshold(tiffImg, removeBackgroundFirst)
% COMPUTE_THRESHOLD Compute the image threshold in each block of the image.
%
% threshold = Compute_Threshold(tiffImg, removeBackgroundFirst)
%
% Input
%   removeBackground : logical flag. If true, and if the TiffImg has a
%     Background, then the background will be removed before computing the
%     image threshold.
%
% Output
%   threshold : Matrix giving the smoothed threshold computed accross the
%     image. Even if no output is requested the threshold is still stored
%     to the tiffImg object.
%
% Note : If tiffImg.Surface_Smoothing_Radius is non-NaN, then the threshold
% will be smoothed with a Lowess smoothing surface.

% James Kapaldo


try
    threshold = zeros(tiffImg.numBlcks,obj.imageClass);
    removeBackgroundFirst = removeBackgroundFirst && ~isempty(tiffImg.BG_fun);
    
    if tiffImg.Verbose
        fprintf('Starting threshold calculation...\n');
    end
    
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose);
    progress.start();
    
    % Iterage over x blocks
    for blck_x = 1:tiffImg.numBlcks(2)
        
        tiffImg.open() % Open Image
        
        % Iterate over y blocks
        for blck_y = 1:tiffImg.numBlcks(1)
            
            % Read in image block
            [I,x,y] = getBlock(tiffImg, blck_x, blck_y);
            
            % Remove background
            if removeBackgroundFirst
                I = I - tiffImg.BG_fun(x,y);
            end
            
            % Filter image
            if tiffImg.Use_GPU
                I = gpuArray(I);
                Is = imfilter(I,tiffImg.Image_Smooth_Kernel,'symmetric'); clear I
                I = gather(Is); clear Is;
            else
                I = imfilter(I,tiffImg.Image_Smooth_Kernel,'symmetric');
            end
            
            % Compute threshold
            threshold(blck_y,blck_x) = otsuthresh_scale(I,'log');
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        progress.iteration_end(); % Update progress counter
    end % x block
    
    % Increase the threshold slightly. The log scale is a bit low.
    threshold = 1.2 * threshold;
    
    % Smooth the threshold surface
    threshold = smoothSurface(tiffImg, threshold);
    
    % Save output for later reference
    tiffImg.threshold = struct('x',tiffImg.xCenters,'y',tiffImg.yCenters,'Z',threshold,'removeBackgroundFirst',removeBackgroundFirst);
    
    % Create function for internally evaluating the threshold in other
    % functions.
    tiffImg.threshold_fun = generateFunction(tiffImg, threshold);
    
    % Modify a flag letting other functions know if they need to apply the
    % threshold after the background correction
    if removeBackgroundFirst
        tiffImg.Threshold_After_Background = true;
    else
        tiffImg.Threshold_After_Background = false;
    end
    
    if tiffImg.Verbose
        fprintf('Threshold calculation finished.\n');
    end
    
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
