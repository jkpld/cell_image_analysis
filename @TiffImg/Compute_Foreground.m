function fg = Compute_Foreground(tiffImg, computeXStripeArtifact)
% COMPUTE_FOREGROUND Compute the image foreground in each block of the
% image.
%
% foreground = Compute_Foreground(tiffImg)
%
% Output
%   foreground : Matrix giving the smoothed foreground computed accross the
%     image. Even if no output is requested the foreground is still stored
%     to the tiffImg object.
%   computeXStripeArtifact : logical flag. If true, then the x stripe
%     artifact will be computed after computing the smoothed foreground.
%     (Default, false)
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

if nargin < 2
    computeXStripeArtifact = false;
end

if isempty(tiffImg.threshold_fun)
    error('Compute_Foreground:noThreshold','The image threshold must be computed before computing the foreground.');
end

try
    hasBackground = ~isempty(tiffImg.BG_smooth);
    if hasBackground
        BG_fun = generateFunction(tiffImg, tiffImg.BG_smooth, tiffImg.BG_Xstripe, true);
    end
    
    FG = zeros(tiffImg.numBlcks, tiffImg.workingClass);
    seD2 = strel('diamond',2);
    
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose, 'name', 'Computing foreground,');
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
                BG = BG_fun(x,y);
                
                if tiffImg.Threshold_After_Background
                    Is = Is - BG;
                    BW = imerode(Is > threshold, seD2);
                else
                    BW = imerode(Is > threshold, seD2);
                    Is = Is - BG;
                end
            else
                BW = imerode(Is > threshold, seD2);
            end
            
            % Get the median background intensity
            FGi = mean(Is(BW),'omitnan'); clear Is BW
            
            if tiffImg.Use_GPU
                FGi = gather(FGi);
            end
            if isnan(FGi)
                FGi = 0;
            end
            
            FG(blck_y,blck_x) = FGi;
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        progress.iteration_end(); % Update progress counter
    end % x block
    
    % Smooth the foreground surface
    FG = smoothSurf(tiffImg, FG);
    
    % Save output for later reference
    tiffImg.FG_smooth = struct('x',tiffImg.xCenters,'y',tiffImg.yCenters,'Z',FG);
    
    if computeXStripeArtifact
        Compute_StripeArtifact(tiffImg,'f');
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