function st = Compute_StripeArtifact(tiffImg, BG_or_FG)
% COMPUTE_STRIPEARTIFACT Compute the X stripe artifact that appears in some
% images due to bad microscope callibration.
%
% stripe = Compute_StripeArtifact(tiffImg)
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
[isBG, str, BG_fun, FG_fun] = parse_input(tiffImg, BG_or_FG);

try
    if tiffImg.Verbose
        fprintf('Starting %s stripe calculation...\n', str);
    end
    
    stripe = zeros(1, tiffImg.imageSize(2), obj.imageClass);
    y = 1:tiffObj.imageSize(1);
    seD2 = strel('diamond',2);
    
    progress = displayProgress(tiffImg.numTiles(2),'number_of_displays', 15,'active',tiffImg.Verbose);
    progress.start();
    
    for row_x = 1:tiffObj.numTiles(2)

        % initialize array for image stripe
        I = zeros(tiffObj.imageSize(1),tiffObj.tileSize(2),tiffImg.imageClass);
        
        % Read in image stripe
        counter = 1;
        tiffObj.open()
        for tile_y = tiffObj.tiles(:,row_x)'
            tmp = tifflib('readEncodedTile',tiffObj.FileID,tile_y);%readEncodedTile(t,cTiles(i));
            [m,n] = size(tmp);
            I(counter:counter+m-1, 1:n) = tmp;
            counter = counter + m;
        end
        tiffObj.close()
        I(:,n+1:end) = [];
        
        % Compute x pixels
        x = (1:size(I,2)) + (row_x-1)*tiffObj.tileSize(2);
                
        % Get threshold for slice
        threshold = tiffImg.threshold_fun(x,y);
        
        if tiffImg.Use_GPU
            I = gpuArray(I);
        end
        
        % Smooth image
        Is = imfilter(I, tiffImg.Image_Smooth_Kernel, 'symmetric'); clear I

        % Remove background before thresholding if supposed to
        if tiffImg.Threshold_After_Background
            Is = Is - BG_fun(x,y);
        end
        
        % Get background/foreground object mask
        % If computing the background stripe than we want to NaN out the
        % foreground, and vise-versa.
        if isBG
            BW = imdilate(Is > threshold, seD2); 
        else
            BW = imdilate(Is < threshold, seD2); 
        end
        
        % Remove background after thresholding if supposed to
        if ~tiffImg.Threshold_After_Background
            Is = Is - BG_fun(x,y);
        end
        
        % Remove foreground if computing foreground stripe
        if ~isBG
            Is = Is - FG_fun(x,y);
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
        tiffImg.BG_stripe = stripe;
    else
        tiffImg.FG_stripe = stripe;
    end
    
    if tiffImg.Verbose
        fprintf('%s stripe calculation finished.\n',str);
    end
    
    % Output the stripe if requested.
    if nargout > 0
        st = stripe;
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


function [BGstripe, str, BG_fun, FG_fun] = parse_input(tiffImg, BG_or_FG)
if isempty(tiffImg.threshold_fun)
    error('Compute_StripeArtifact:noThreshold','The image threshold and then the background or foreground must be computed before computing the stripe artifact.');
end

idx = find(strncmpi(BG_or_FG,{'background','foreground'},numel(BG_or_FG)));
if numel(idx) > 1
    error('Compute_StripeArtifact:ambigiousInput','The BG_or_FG string input must be an unambigious match to either "background" or "foreground".');
else
    BGstripe = idx == 0;
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
    % the smooth background.
    BG_fun = generateFunction(tiffImg,tiffImg.BG_smooth);
    FG_fun = [];
else
    str = 'foreground';
    % If computing the foreground stripe, then we need a function to
    % evaluate the full background (if any)
    if ~isempty(tiffImg.BG_smooth)
        BG_fun = tiffImg.BackgroundEvalFun;
    else
        BG_fun = @(~,~) 0;
    end
    
    % We also need a function to evaluate the smoothe foreground.
    FG_fun = generateFunction(tiffImg,tiffImg.FG_smooth);
end    
end