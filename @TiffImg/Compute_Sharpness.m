function sh = Compute_Sharpness(tiffImg)
% COMPUTE_SHARPNESS Compute the image sharpness in each block of the
% image.
%
% sharpness = Compute_Sharpness(tiffImg)
%
% Output
%   sharpness : Matrix giving the smoothed sharpness computed accross the
%     image. Even if no output is requested the sharpness is still stored
%     to the tiffImg object.
%
% Note : The sharpness is computed by taking the fourier transform and
% counting the number of values above some threshold.
%
% Note : If tiffImg.Surface_Smoothing_Radius is non-NaN, then the sharpness
% will be smoothed with a Lowess smoothing surface.

% James Kapaldo

try   
    sharp = zeros(tiffImg.numBlcks, tiffImg.workingClass);
    
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose, 'name', 'Computing image sharpness,');
    progress.start();
    
    % Iterate over x blocks
    for blck_x = 1:tiffImg.numBlcks(2)
        
        tiffImg.open() % Open image
        
        % Iterate over y blocks
        for blck_y = 1:tiffImg.numBlcks(1)
            
            % Get image block
            I = getBlock(tiffImg,blck_x,blck_y);
            N = numel(I); % number of pixels
            
            if tiffImg.Use_GPU
                I = gpuArray(I);
            end
            
            % Smooth image
            Is = imfilter(I, tiffImg.Image_Smooth_Kernel, 'symmetric'); clear I
            
            % Compute fft
            IsFFT = abs(fft2(Is)); clear Is
            IsFFT_flat = IsFFT(:); clear IsFFT
            
            % Compute sharpness
            sharp_i = sum(IsFFT_flat > max(IsFFT_flat)/1000); clear IsFFT_flat
            
            if tiffImg.Use_GPU
                sharp_i = gather(sharp_i);
            end
            
            sharp(blck_y,blck_x) = sharp_i / sqrt(N); % Normalize sharpness by sqrt of number of pixels.
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        progress.iteration_end(); % Update progress counter
    end % x block
    
    % Smooth the sharpness surface
    sharp = smoothSurf(tiffImg, sharp);
    
    % Save output for later reference
    tiffImg.Sharpness = struct('x',tiffImg.xCenters,'y',tiffImg.yCenters,'Z',sharp);
    
    % Output the sharpness matrix if requested.
    if nargout > 0
        sh = sharp;
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