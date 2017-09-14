function fg = Compute_Foreground(tiffImg)
% COMPUTE_FOREGROUND Compute the image foreground in each block of the
% image.
%
% foreground = Compute_Foreground(tiffImg)
%
% Output
%   foreground : Matrix giving the smoothed foreground computed accross the
%     image. Even if no output is requested the foreground is still stored
%     to the tiffImg object.
%
% Note : Computing the foreground requires that the image threshold has
% already been calculated. The foreground is calculated by first
% thresholding the image to obtain the foreground regions, and then
% computing the median of the foreground region. If a background has
% already been computed, then the background is taken into account before
% computing the foreground.
%
% Note : If tiffImg.Surface_Smoothing_Radius is non-NaN, then the
% foreground will be smoothed with a Lowess smoothing surface.

% James Kapaldo

if isempty(tiffImg.threshold_fun)
    error('Compute_Foreground:noThreshold','The image threshold must be computed before computing the foreground.');
end

try
    hasBackground = ~isempty(tiffImg.BG_fun);
    
    if tiffImg.Verbose
        fprintf('Starting foreground calculation...\n');
    end
    
    FG = zeros(tiffImg.numBlcks, obj.imageClass);
    seD2 = strel('diamond',2);
    
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose);
    progress.start();
    
    % Iterate over x blocks
    for blck_x = 1:tiffImg.numBlcks(2)
        
        tiffImg.open() % Open image
        
        % Iterate over y blocks
        for blck_y = 1:tiffImg.numBlcks(1)
            
            % Get image block
            [I,x,y] = getBlock(tiffImg,blck_x,blck_y);
            
            % Get threshold for block.
            threshold = tiffImg.threshold_fun(x,y);

            if tiffImg.Use_GPU
                I = gpuArray(I);
            end
            
            % Smooth image
            Is = imfilter(I, tiffImg.Image_Smooth_Kernel, 'symmetric'); clear I
            
            % Get the foreground mask
            % - erode the foreground mask so that it is farther away from
            % the background
            if hasBackground
                BG = tiffImg.BG_fun(x,y);
                
                % The below code is long, but it reduces the amount of
                % memory used by the GPU.
                if tiffImg.Threshold_After_Background
                    if tiffImg.Use_GPU
                        BG = gpuArray(BG);
                    end
                    Is = Is - BG; clear BG
                    if tiffImg.Use_GPU
                        threshold = gpuArray(threshold);
                    end
                    BW = imerode(Is > threshold, seD2); clear threshold
                else
                    if tiffImg.Use_GPU
                        threshold = gpuArray(threshold);
                    end
                    BW = imerode(Is > threshold, seD2); clear threshold

                    if tiffImg.Use_GPU
                        BG = gpuArray(BG);
                    end
                    Is = Is - BG; clear BG
                end
            else
                if tiffImg.Use_GPU
                    threshold = gpuArray(threshold);
                end
                BW = imerode(Is > threshold, seD2); clear threshold
            end
            
            % Get the median background intensity
            FGi = median(Is(BW)); clear Is BW
            
            if tiffImg.Use_GPU
                FGi = gather(FGi);
            end
            
            FG(blck_y,blck_x) = FGi;
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        progress.iteration_end(); % Update progress counter
    end % x block
    
    % Smooth the foreground surface
    FG = smoothSurface(tiffImg, FG);
    meanFG = mean(FG(:));
    
    % Save output for later reference
    tiffImg.Foreground_smooth = struct('x',tiffImg.xCenters,'y',tiffImg.yCenters,'Z',FG);
    tiffImg.Foreground_mean = meanFG;
    
    % Create function for internally evaluating the foreground in other
    % functions.
    tiffImg.FG_smooth_fun = generateFunction(tiffImg, FG);
    tiffImg.FG_mean = meanFG;
    tiffImg.FG_fun = generateFunction(tiffImg, FG, [], -meanFG);
    
    if tiffImg.Verbose
        fprintf('Foreground calculation finished.\n');
    end
    
    % Output the background matrix if requested.
    if nargout > 0
        fg = FG;
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