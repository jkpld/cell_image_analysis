function st = Compute_StripeArtifact(tiffImg, BG_or_FG, varargin)
% COMPUTE_STRIPEARTIFACT Compute the X stripe artifact that appears in some
% images due to bad microscope callibration.
%
% stripe = Compute_StripeArtifact(tiffImg)
%
% Input
%   BG_or_FG : String that unambigously matches either 'background' or
%     'foreground'. This will determine if the background stripe or the
%     foreground stripe is computed.
%
% Optional param/value Input
%   'Object_Mask' : a TiffImg object that gives the object mask. If
%     provided, then the threshold will not need to be computed.
%
% Output
%   stripe : 1 x nCols array giving the stripe artifact computed accross
%     the y direction. nCols is the number of image columns. Even if no
%     output is requested the stripe is still stored to the tiffImg object.
%
% Note : Computing the stripe artifact requires that the image threshold
% has already been calculated and that the smooth background or smooth
% foreground has been calculated (depending on if the stripe is to be
% calculated for the background or foreground). The smooth back/fore-ground
% is removed and then the median pixel value accross entire columns is
% computed; these values make the final stripe artifact.
%
% Note : If computing the foreground stripe, then any background that was
% previously computed is removed before computing the foreground.

% James Kapaldo

% Check input and get required functions
[isBG, str, BG_fun, FG_fun, object_mask] = parse_input(tiffImg, BG_or_FG, varargin{:});
Use_Mask = ~isEmpty(object_mask);
if Use_Mask, object_mask.blockSize = tiffImg.blockSize; end

try    
    stripe = zeros(1, tiffImg.imageSize(2), tiffImg.workingClass);
    y = (1:tiffImg.imageSize(1)).';
    seD2 = strel('diamond',2);
    
    progress = displayProgress(tiffImg.numTiles(2),'number_of_displays', 15,'active',tiffImg.Verbose, 'name',['Computing ', str, ' stripe,']);
    progress.start();
    
    for col_x = 1:tiffImg.numTiles(2)

        % Get column
        [I,x] = getColumn(tiffImg, col_x);
        
        if tiffImg.Use_GPU
            I = gpuArray(I);
        end
        
        % Smooth image
        Is = imfilter(I, tiffImg.Image_Smooth_Kernel, 'symmetric'); clear I
        
        if Use_Mask
            BWo = getColumn(object_mask, col_x);
                
            if tiffImg.Use_GPU
                BWo = gpuArray(BWo);
            end
            
            if ~isBG
                BWo = ~BWo;
            end

            BW = imdilate(BWo, seD2); clear BWo
        else
            
            % Get threshold for slice
            threshold = tiffImg.threshold_fun(x,y);

            % Remove background before thresholding if supposed to
            if tiffImg.Threshold_After_Background
                Ist = Is - BG_fun(x,y);
            else
                Ist = Is;
            end

            % Get background/foreground object mask
            % If computing the background stripe than we want to NaN out the
            % foreground, and vise-versa.
            if isBG
                BW = imdilate(Ist > threshold, seD2); 
            else
                BW = imdilate(Ist < threshold, seD2); 
            end
        end
        
        % Divide the image by the background or foreground. This will make
        % the x-stripe a multiplicitive factor.
        if isBG
            Is = Is ./ BG_fun(x,y);
        else
            if ~tiffImg.Threshold_After_Background
                Is = Is - BG_fun(x,y);
            end

            Is = Is ./ FG_fun(x,y);
        end
        
        Is(BW) = nan; clear BW
        stripe_i = median(Is, 1, 'omitnan');
        
        if tiffImg.Use_GPU
            stripe_i = gather(stripe_i);
        end
        stripe(x) = stripe_i;
        
        progress.iteration_end(); % Update progress counter
    end

    
    % Save stripe
    if isBG
        tiffImg.BG_Xstripe = stripe;
    else
        tiffImg.FG_Xstripe = stripe;
    end
    
    % Output the stripe if requested.
    if nargout > 0
        st = stripe;
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


function [BGstripe, str, BG_fun, FG_fun, object_mask] = parse_input(tiffImg, BG_or_FG, varargin)

p = inputParser;
p.FunctionName = 'Compute_StripeArtifact';

addParameter(p,'Object_Mask', TiffImg(), @(t) isa(t,'TiffImg'));

parse(p,varargin{:})
object_mask = p.Results.Object_Mask;  

if isempty(tiffImg.threshold_fun) && isempty(object_mask)
    error('Compute_StripeArtifact:noThreshold','There must be an Object_Mask, or the image threshold and then the background or foreground must be computed before computing the stripe artifact.');
end

idx = find(strncmpi(BG_or_FG,{'background','foreground'},numel(BG_or_FG)));

if numel(idx) > 1
    error('Compute_StripeArtifact:ambigiousInput','The BG_or_FG string input must be an unambigious match to either "background" or "foreground".');
else
    BGstripe = idx == 1;
end

if BGstripe
    if isempty(tiffImg.BG_smooth)
        error('Compute_StripeArtifact:missingBackground','The background must be computed before computing the background stripe artifact.');
    end
else
    if isempty(tiffImg.FG_smooth)
        error('Compute_StripeArtifact:missingForeground','The foreground must be computed before computing the foreground stripe artifact.');
    end
end

if BGstripe
    str = 'background';
    % If computing bacgkround stripe, then we need a function to evaluate
    % the smooth background. This will automatically offset the surface by
    % the surface mean.
    BG_fun = generateFunction(tiffImg,tiffImg.BG_smooth);
    FG_fun = [];
else
    str = 'foreground';
    % If computing the foreground stripe, then we need a function to
    % evaluate the full background (if any). This will automatically offset
    % the surface by the surface mean.
    if ~isempty(tiffImg.BG_smooth)
        BG_fun = generateFunction(tiffImg,tiffImg.BG_smooth,tiffImg.BG_Xstripe,true);
    else
        BG_fun = @(~,~) 0;
    end
    
    % We also need a function to evaluate the smoothe foreground.
    FG_fun = generateFunction(tiffImg,tiffImg.FG_smooth);
end




end
