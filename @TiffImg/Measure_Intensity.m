function feature = Measure_Intensity(tiffImg, Use_Parallel, varargin)
% MEASURE_INTENSITY Compute the integrated intensity for each object in the
% given image.
%
% feature = Measure_BasicProps(tiffImg, Use_Parallel)
%
% Input
%   Use_Parallel : logical flag.
%
% Optional param/value Input
%   'Object_Mask' : a TiffImg object that gives the object mask. If
%     provided, then the threshold will not need to be computed.
%
% Output
%   feature : N x 1 array, where N is the number of objects, giving the
%     integrated intensity.

% James Kapaldo

Use_Parallel = logical(Use_Parallel);
object_mask = parse_input(varargin{:});
Use_Mask = ~isEmpty(object_mask);
if Use_Mask, object_mask.blockSize = tiffImg.blockSize; end

if isempty(tiffImg.threshold_fun) && ~Use_Mask
    error('Compute_Foreground:noThreshold','The image threshold must be computed before computing the basic object properties, or an Object_Mask must be provided.');
end

try
    if ~tiffImg.Threshold_After_Correction
        CorrectionFunction = generateCorrectionFunction(tiffImg);
    end
    
    % Initialize features
    feature = zeros(50000, 1, 'single');
    feature_count = 1;
    
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose, 'name', 'Computing object integrated intensities,');
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
                BW = getBlock(object_mask, blck_x, blck_y);
                
                if tiffImg.Use_GPU
                    BW = gpuArray(BW);
                end
                
                % apply image corrections
                Is = CorrectionFunction(Is,x,y);
            else
                % Get threshold for block.
%                 if ~tiffImg.Threshold_After_Correction
                    threshold = tiffImg.threshold_fun(x,y);
%                 end

                % Apply corrections if before threshold
                if tiffImg.Threshold_After_Correction
                    Is = tiffImg.Threshold_Correction(Is,x,y);
                    BW = Is > threshold;
                    if ~isempty(tiffImg.Secondary_Correction)
                        Is = tiffImg.Secondary_Correction(Is,x,y);
                    end
                else
                    BW = Is > threshold;
                    Is = CorrectionFunction(Is,x,y);
                end
            end
            
            % Label matrix            
            BW = clear2dborder(BW);
            
            % If using gpu, first gather them before sending them to
            % bp.Compute()
            
            if tiffImg.Use_GPU
                I = gather(Is); clear Is
                BWc = gather(BW); clear BW
            else
                I = Is; clear Is
                BWc = BW; clear BW
            end
            
            
            if ~any(BWc(:))
                continue;
            else
                CC = bwconncomp(BWc);
                N = CC.NumObjects;
                pixList = CC.PixelIdxList;
                x_i = zeros(N,1,'single');
                if Use_Parallel
                    I = cellfun(@(x) I(x), pixList, 'UniformOutput',false);
                    parfor i = 1:N
                        x_i(i) = sum(I{i});
                    end
                else
                    for i = 1:N
                        x_i(i) = sum(I(pixList{i}));
                    end
                end                

                if N+feature_count-1 > size(feature,1)
                    feature = [feature; zeros(50000, 1,'single')]; %#ok<AGROW>
                end
                feature(feature_count:N+feature_count-1,:) = x_i;
                feature_count = feature_count + N;
            end
            
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        if Use_Mask, object_mask.close(); end
        progress.iteration_end(); % Update progress counter
    end % x block
    
    % Remove extra allocated features
    feature(feature_count:end,:) = [];
    
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
p.FunctionName = 'Measure_Intensity';

addParameter(p,'Object_Mask', TiffImg(), @(t) isa(t,'TiffImg'));

parse(p,varargin{:})
object_mask = p.Results.Object_Mask;        
end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
