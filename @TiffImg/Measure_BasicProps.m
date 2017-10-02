function [feature,feature_names] = Measure_BasicProps(tiffImg, channelName, varargin)
% MEASURE_BASICPROPS Compute the centroid, area, and integrated intensity
% for each object in the given image.
%
% [feature,feature_names] = Measure_BasicProps(tiffImg, channelName)
%
% Input
%   channelName : The name of the input channel. This is only used to
%     produce proper feature_Names.
%
% Optional param/value Input
%   'Object_Mask' : a TiffImg object that gives the object mask. If
%     provided, then the threshold will not need to be computed.
%
% Output
%   feature : N x 4 array where N is the number of objects and 4 is the
%     number of features (Centroid_x, Centroid_y, Area, Intensity)
%   feature_names : 1 x 4 string array giving the name of each feature

% James Kapaldo

object_mask = parse_input(varargin{:});
Use_Mask = ~isEmpty(object_mask);
if Use_Mask, object_mask.blockSize = tiffImg.blockSize; end

if nargin < 2
    channelName = '';
end
if isempty(tiffImg.threshold_fun) && ~Use_Mask
    error('Compute_Foreground:noThreshold','The image threshold must be computed before computing the basic object properties, or an Object_Mask must be provided.');
end

try
    CorrectionFunction = generateCorrectionFunction(tiffImg);
    ThresholdCorrection_isCurrent = isequal(tiffImg.Threshold_CorrectionsExpression, tiffImg.Current_Image_Correction_Expression);
    
    if tiffImg.Threshold_After_Correction && ~ThresholdCorrection_isCurrent
        Threshold_Correction = generateFunctionFromExpression(tiffImg, tiffImg.Threshold_CorrectionsExpression, true);
    end
    
    if tiffImg.Threshold_After_Correction && ~Use_Mask
        threshold = median(tiffImg.threshold.Z(:));
    end
    
    % Basic props measurer
    bp = BasicProps(channelName,struct('Use_GPU',tiffImg.Use_GPU));
    
    % Initialize features
    feature_names = bp.FeatureNames;
    num_feature = numel(feature_names);
    feature = zeros(50000, num_feature, 'single');
    feature_count = 1;
    
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose, 'name', 'Computing basic properties,');
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
                if ~tiffImg.Threshold_After_Correction
                    threshold = tiffImg.threshold_fun(x,y);
                end

                % Apply corrections if before threshold
                if tiffImg.Threshold_After_Correction
                    if ThresholdCorrection_isCurrent
                        Is = CorrectionFunction(Is,x,y);
                        BW = Is > threshold;
                    else
                        BW = Threshold_Correction(Is,x,y) > threshold;
                        Is = CorrectionFunction(Is,x,y);
                    end
                else
                    BW = Is > threshold;
                    Is = CorrectionFunction(Is,x,y);
                end
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
p.FunctionName = 'Measure_BasicProps';

addParameter(p,'Object_Mask', TiffImg(), @(t) isa(t,'TiffImg'));

parse(p,varargin{:})
object_mask = p.Results.Object_Mask;        
end