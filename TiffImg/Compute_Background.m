function bg = Compute_Background(tiffImg)
% COMPUTE_BACKGROUND Compute the image background in each block of the
% image.
%
% background = Compute_Background(tiffImg)
%
% Output
%   background : Matrix giving the smoothed background computed accross the
%     image. Even if no output is requested the background is still stored
%     to the tiffImg object.
%
% Note : Computing the background requires that the image threshold has
% already been calculated. The background is calculated by first
% thresholding the image to obtain the background regions, and then
% computing the median of the background region.
%
% Note : If tiffImg.Surface_Smoothing_Radius is non-NaN, then the
% background will be smoothed with a Lowess smoothing surface.

% James Kapaldo

if isempty(tiffImg.threshold_fun)
    error('Compute_Background:noThreshold','The image threshold must be computed before computing the background.');
end

try
    if tiffImg.Verbose
        fprintf('Starting background calculation...\n');
    end
    
    BG = zeros(tiffImg.numBlcks, obj.imageClass);
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
                threshold = gpuArray(threshold);
            end
            
            % Smooth image
            Is = imfilter(I, tiffImg.Image_Smooth_Kernel, 'symmetric'); clear I
            
            % Get the background mask
            % - erode the background mask so that it is farther away from
            % the objects
            BW = imerode(Is < threshold, seD2); clear threshold
            
            % Get the median background intensity
            BGi = median(Is(BW)); clear Is BW
            
            if tiffImg.Use_GPU
                BGi = gather(BGi);
            end
            
            BG(blck_y,blck_x) = BGi;
            
            %             Th = thrsh_fun({y,xt});
            %             Th_d = gpuArray(Th);
            %
            %             I = single(I);
            %             I_d = gpuArray(I);
            %             Is_d = imfilter(I_d,H,'symmetric');
            %
            %             % Background and foreground mask
            %             BW_d = Is_d > Th_d;
            %             clear Th_d
            %             BW_BG_d = ~imdilate(BW_d,seD2);
            %             BW_FG_d = imerode(BW_d,seD2);
            %
            %             % Background and forground value
            %             BGt_d = median(Is_d(BW_BG_d));
            %             FGt_d = median(Is_d(BW_FG_d));
            %             BG(blck_y,blck_x) = gather(BGt_d);
            %             FG(blck_y,blck_x) = gather(FGt_d);
            %             clear I_d BGt_d FGt_d BW_BG_d BW_FG_d BW_d
            %
            %             % Bluriness
            %             N = numel(I);
            %             IsFFT_d = abs(fft2(Is_d));
            %             BLUR(blck_y,blck_x) = gather(sum(IsFFT_d(:) > max(IsFFT_d(:))/1000)) / sqrt(N);
            %             clear Is_d IsFFT_d
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        progress.iteration_end(); % Update progress counter
    end % x block
    
    % Smooth the background surface
    BG = smoothSurface(tiffImg, BG);
    meanBG = mean(BG(:));
    
    % Save output for later reference
    tiffImg.Background_smooth = struct('x',tiffImg.xCenters,'y',tiffImg.yCenters,'Z',BG);
    tiffImg.Background_mean = meanBG;
    
    % Create function for internally evaluating the background in other
    % functions.
    tiffImg.BG_smooth_fun = generateFunction(tiffImg, BG);
    tiffImg.BG_mean = meanBG;
    tiffImg.BG_fun = generateFunction(tiffImg, BG, [], -meanBG);
    
    if tiffImg.Verbose
        fprintf('Background calculation finished.\n');
    end
    
    % Output the background matrix if requested.
    if nargout > 0
        bg = BG;
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