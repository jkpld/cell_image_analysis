function bg = Compute_Background(tiffImg, computeXStripeArtifact, varargin)
% COMPUTE_BACKGROUND Compute the image background in each block of the
% image.
%
% background = Compute_Background(tiffImg)
%
% Input
%   computeXStripeArtifact : (optional) logical flag. If true, then the x
%     stripe artifact will be computed after computing the smoothed
%     background. (Default false)
%
% Optional param/value Input
%   'Object_Mask' : a TiffImg object that gives the object mask. If
%     provided, then the threshold will not need to be computed.
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

object_mask = parse_input(varargin{:});
Use_Mask = ~isEmpty(object_mask);
if Use_Mask, object_mask.blockSize = tiffImg.blockSize; end

if nargin < 2
    computeXStripeArtifact = false;
end

if isempty(tiffImg.threshold_fun) && ~Use_Mask
    error('Compute_Background:noThreshold','The image threshold must be computed before computing the background, or an Object_Mask must be provided.');
end

try    

    BG = zeros(tiffImg.numBlcks, tiffImg.workingClass);
    seD2 = strel('diamond',2);
    
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose, 'name', 'Computing background,');
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
                
                BW = imerode(~BWo, seD2); clear BWo
            else
                % Get threshold for block.
                threshold = tiffImg.threshold_fun(x,y);

                if tiffImg.Use_GPU
                    threshold = gpuArray(threshold);
                end               

                % Get the background mask
                % - erode the background mask so that it is farther away from
                % the objects
                BW = imerode(Is < threshold, seD2); clear threshold
            end
            
            % Get the median background intensity
            BGi = median(Is(BW),'omitnan'); clear Is BW
            
            if tiffImg.Use_GPU
                BGi = gather(BGi);
            end
            if isnan(BGi)
                BGi = 0;
            end
            BG(blck_y,blck_x) = BGi;
            
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        if Use_Mask, object_mask.close(); end
        progress.iteration_end(); % Update progress counter
    end % x block
    
    % Smooth the background surface
    BG = smoothSurf(tiffImg, BG);
    
    % Save background
    tiffImg.BG_smooth = struct('x',tiffImg.xCenters,'y',tiffImg.yCenters,'Z',BG);
    
    if computeXStripeArtifact
        Compute_StripeArtifact(tiffImg,'b', 'Object_Mask', object_mask);
    end
    
    % Output the background matrix if requested.
    if nargout > 0
        bg = BG;
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
p.FunctionName = 'Compute_Background';

addParameter(p,'Object_Mask', TiffImg(), @(t) isa(t,'TiffImg'));

parse(p,varargin{:})
object_mask = p.Results.Object_Mask;        
end