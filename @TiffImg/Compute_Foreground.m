function fg = Compute_Foreground(tiffImg, computeXStripeArtifact, varargin)
% COMPUTE_FOREGROUND Compute the image foreground in each block of the
% image.
%
% foreground = Compute_Foreground(tiffImg)
%
% Input
%   computeXStripeArtifact : logical flag. If true, then the x stripe
%     artifact will be computed after computing the smoothed foreground.
%     (Default, false)
%
% Optional param/value Input
%   'Object_Mask' : a TiffImg object that gives the object mask. If
%     provided, then the threshold will not need to be computed.
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

object_mask = parse_input(varargin{:});
Use_Mask = ~isEmpty(object_mask);
if Use_Mask, object_mask.blockSize = tiffImg.blockSize; end

if nargin < 2
    computeXStripeArtifact = false;
end

if isempty(tiffImg.threshold_fun) && ~Use_Mask
    error('Compute_Foreground:noThreshold','The image threshold must be computed before computing the foreground, or an Object_Mask must be provided.');
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
        if Use_Mask, object_mask.open(); end
        
        % Iterate over y blocks
        for blck_y = 1:tiffImg.numBlcks(1)
            
            % Get image block
            [I,x,y] = getBlock(tiffImg,blck_x,blck_y);
            
            if tiffImg.Use_GPU
                I = gpuArray(I);
            end
            
            % Smooth image
            Is = imfilter(I, tiffImg.Image_Smooth_Kernel, 'symmetric'); clear I
            
            if Use_Mask
                BWo = getBlock(object_mask, blck_x, blck_y);
                
                if tiffImg.Use_GPU
                    BWo = gpuArray(BWo);
                end
                
                BW = imerode(BWo, seD2); clear BWo
                
                if hasBackground
                    BG = BG_fun(x,y);
                    Is = Is - BG; clear BG
                end
            else
                
                % Get threshold for block.
                threshold = tiffImg.threshold_fun(x,y);

                % Get the foreground mask
                % - erode the foreground mask so that it is farther away from
                % the background
                if hasBackground
                    BG = BG_fun(x,y);

                    if tiffImg.Threshold_After_Background
                        Is = Is - BG; clear BG
                        BW = imerode(Is > threshold, seD2);
                    else
                        BW = imerode(Is > threshold, seD2);
                        Is = Is - BG; clear BG
                    end
                else
                    BW = imerode(Is > threshold, seD2);
                end
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
        if Use_Mask, object_mask.close(); end
        progress.iteration_end(); % Update progress counter
    end % x block
    
    % Smooth the foreground surface
    FG = smoothSurf(tiffImg, FG);
    
    % Save output for later reference
    tiffImg.FG_smooth = struct('x',tiffImg.xCenters,'y',tiffImg.yCenters,'Z',FG);
    
    if computeXStripeArtifact
        Compute_StripeArtifact(tiffImg,'f', 'Object_Mask', object_mask);
    end
    
    % Output the background matrix if requested.
    if nargout > 0
        fg = FG;
    end
    
    
catch ME
    tiffImg.close();
    if Use_Mask, object_mask.close(); end
    if tiffImg.Use_GPU
        gpuDevice([]);
    end
    fprintf('entered catch statement\n')
    rethrow(ME)
end

tiffImg.close();
if Use_Mask, object_mask.close(); end

end % function

function object_mask = parse_input(varargin)
p = inputParser;
p.FunctionName = 'Compute_Foreground';

addParameter(p,'Object_Mask', TiffImg(), @(t) isa(t,'TiffImg'));

parse(p,varargin{:})
object_mask = p.Results.Object_Mask;        
end