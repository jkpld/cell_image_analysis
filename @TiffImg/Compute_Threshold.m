function th = Compute_Threshold(tiffImg, applyCorrectionsFirst)
% COMPUTE_THRESHOLD Compute the image threshold in each block of the image.
%
% threshold = Compute_Threshold(tiffImg, removeBackgroundFirst, removeForegroundFirst)
%
% Input
%   applyCorrectionsFirst : logical flag. If true, and then the image will
%     be corrected using the Current_Image_Correction_Expression before
%     computing the threshold
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
    applyCorrectionsFirst = false;
end

try
    threshold = zeros(tiffImg.numBlcks,tiffImg.workingClass);
    
    if applyCorrectionsFirst
        CorrectionFunction = generateCorrectionFunction(tiffImg);
        % Offset by the median background value
%         if ~isempty(tiffImg.BG_offset)
%             mBG = median(tiffImg.BG_offset.Z(:))
%             CorrectionFunction = @(I,x,y) CorrectionFunction(I,x,y) + mBG;
%         end
        scale = 'linear';
    else
        scale = 'log';
    end
%         scale = 'log';
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose,'name',sprintf('Computing threshold (%s),',scale));
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
            if applyCorrectionsFirst
                I = CorrectionFunction(I,x,y);
            end

            % Compute threshold
%             if any(I(:)>1)
%                 fprintf('larger than 1!\n')
%             end
% figure
% histogram(I)
            threshold(blck_y,blck_x) = tiffImg.otsuthresh_scale(I,scale);
%             threshold(blck_y,blck_x)
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        progress.iteration_end(); % Update progress counter
    end % x block
    
%     error(' some error')
    
    % Increase the threshold slightly. The log scale is a bit low.
    if strcmp(scale,'log')
        threshold = 1.2 * threshold;
    end    
    threshold = medfilt2(threshold,[3,3],'symmetric');
    % Smooth the threshold surface
    threshold = smoothSurf(tiffImg, threshold);
    
    % Save output for later reference
    tiffImg.threshold = struct('x',tiffImg.xCenters,'y',tiffImg.yCenters,'Z',threshold,'applyCorrectionsFirst',applyCorrectionsFirst);
    
    % Create function for internally evaluating the threshold in other
    % functions.
    tiffImg.threshold_fun = interpolator2d(tiffImg.threshold.x,tiffImg.threshold.y,tiffImg.threshold.Z);
    
    % Modify a flag letting other functions know if they need to apply the
    % threshold after the image correction.
    tiffImg.Threshold_After_Correction = applyCorrectionsFirst;
    % Save the current correction expression
    tiffImg.Threshold_CorrectionsExpression = tiffImg.Current_Image_Correction_Expression;
    if applyCorrectionsFirst
        tiffImg.Threshold_Correction = generateCorrectionFunction(tiffImg);
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
