function [feature,feature_names] = Measure_BasicProps(tiffImg, channelName)
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
    channelName = 'DAPI';
end
if isempty(tiffImg.threshold_fun)
    error('Compute_Foreground:noThreshold','The image threshold must be computed before computing the basic object properties.');
end

try
    hasBackground = ~isempty(tiffImg.BG_smooth);
    if hasBackground
        BG_fun = generateFunction(tiffImg, tiffImg.BG_smooth, tiffImg.BG_Xstripe, true);
    end
    hasForeground = ~isempty(tiffImg.FG_Xstripe);
    if hasForeground       
        FG_fun = generateFunction(tiffImg, tiffImg.FG_smooth, tiffImg.FG_Xstripe, true);
    end
    
    backgroundFirst = hasBackground && tiffImg.Threshold_After_Background;
    backgroundAfter = hasBackground && ~tiffImg.Threshold_After_Background;
    foregroundFirst = hasForeground && tiffImg.Threshold_After_Background;
    foregroundAfter = hasForeground && ~tiffImg.Threshold_After_Background;
    
    % Basic props measurer
    bp = BasicProps(channelName,struct('Use_GPU',tiffImg.Use_GPU));
    
    % Initialize features
    feature_names = bp.FeatureNames;
    num_feature = numel(feature_names);
    feature = zeros(100000, num_feature, 'single');
    feature_count = 1;
    
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose, 'name', 'Computing basic properties,');
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
            
%             if tiffImg.Use_GPU
%                 I = gather(Is); clear Is
%             else
%                 I = Is; clear Is
%             end
            
            % Get the foreground mask
            % - erode the foreground mask so that it is farther away from
            % the background
            if backgroundFirst
                Is = Is - BG_fun(x,y);
            end
            if foregroundFirst
                Is = Is ./ FG_fun(x,y);
            end
            
            BW = Is > threshold;
            
            if backgroundAfter
                Is = Is - BG_fun(x,y);
            end
            if foregroundAfter
                Is = Is ./ FG_fun(x,y);
            end
            
            % Label matrix            
            BW = clear2dborder(BW);
            Ld = bwlabel(BW); clear BW
            
            % If using gpu, first gather them before sending them to
            % bp.Compute()
            
            if tiffImg.Use_GPU
                I = gather(Is); clear Is
                L = gather(Ld); clear Ld
            else
                I = Is; clear Is
                L = Ld; clear Ld
            end
            
            if all(L(:)==0)
                continue;
            else
                x_i = bp.Compute(I, L);
                x_i(:,1) = x_i(:,1) + (blck_x-1)*tiffImg.blockSize;
                x_i(:,2) = x_i(:,2) + (blck_y-1)*tiffImg.blockSize;
                N = size(x_i,1);

                if N+feature_count-1 > size(feature,1)
                    feature = [feature; zeros(50000, num_feature,'single')]; %#ok<AGROW>
                end
                feature(feature_count:N+feature_count-1,:) = x_i;
                feature_count = feature_count + N;
            end
            
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        progress.iteration_end(); % Update progress counter
    end % x block
    
    % Remove extra allocated features
    feature(feature_count:end,:) = [];
    
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