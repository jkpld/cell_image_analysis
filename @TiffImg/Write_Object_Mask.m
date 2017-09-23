function Write_Object_Mask(tiffImg, fileName, partitioner)
% WRITE_OBJECT_MASK Write the (optionally partitioned) object mask to a
% file
%
% Write_Object_Mask(tiffImg, fileName, partitioner)
%
% Input
%   fileName : Name of the file where the mask will be written
%   partitioner : An ObjectPartitioner. If not provided, then the objects
%     will not be partitioned.

% James Kapaldo

if nargin < 3
    partitioner = None;
else
    if ~isa(partitioner, 'ObjectPartitioner')
        error('Write_Object_Mask:badPartitioner', 'The partitioner object must be of class ObjectPartitioner.');
    end
end

% Fields to copy
paramNames = {'SubFileType', ...
    'Photometric', ...
    'PlanarConfiguration', ...
    'ImageLength', ...
    'ImageWidth', ...
    'TileLength', ...
    'TileWidth', ...
    'BitsPerSample', ...
    'SamplesPerPixel', ...
    'SampleFormat', ...
    'Compression', ...
    'ImageDescription', ...
    'Orientation'};

maskID = tifflib('open',char(fileName),'w');

try
    
    % Initialize the tiff tags for the mask image
    tiffImg.open();
    t_ID = tiffImg.FileID;

    for i = 1:numel(paramNames)
        tifflib('setField',maskID,Tiff.TagID.(char(paramNames{i})), ...
            tifflib('getField',t_ID,Tiff.TagID.(char(paramNames{i}))));
    end

    tifflib('setField',maskID,Tiff.TagID.Photometric,4)
    tifflib('setField',maskID,Tiff.TagID.BitsPerSample,1) % set to logical

    tifflib('setField',maskID,Tiff.TagID.ImageDescription, ...
            [tifflib('getField',t_ID,Tiff.TagID.ImageDescription),'|imageMask']);
    
    
    hasBackground = ~isempty(tiffImg.BG_smooth);
    if hasBackground
        BG_fun = generateFunction(tiffImg, tiffImg.BG_smooth, tiffImg.BG_Xstripe, true);
    end
    hasForeground = ~isempty(tiffImg.FG_Xstripe);
    if hasForeground       
        FG_fun = generateFunction(tiffImg, tiffImg.FG_smooth, tiffImg.FG_Xstripe, true);
    end

    backgroundFirst = hasBackground && tiffImg.Threshold_After_Background;
    foregroundFirst = hasForeground && tiffImg.Threshold_After_Background;
    
    % Only need to do corrections afterwards if the partitioner needs the
    % image
    partitionerNeedsImage = ~isa(partitioner,'None') && ...
        partitioner.Restricted_Partitioning && ...
        (~isempty(partitioner.Area_Normalizer) || ...
        ~isempty(partitioner.Intensity_Normalizer));
    backgroundAfter = hasBackground && ~tiffImg.Threshold_After_Background && partitionerNeedsImage;
    foregroundAfter = hasForeground && ~tiffImg.Threshold_After_Background && partitionerNeedsImage;

    if backgroundFirst && foregroundFirst
        threshold = median(tiffImg.threshold.Z(:));
    end
    
    progress = displayProgress(tiffImg.numBlcks(2),'number_of_displays', 15,'active',tiffImg.Verbose, 'name', 'Writing object mask,');
    progress.start();
    
    % Iterate over x blocks
    for blck_x = 1:tiffImg.numBlcks(2)
        
        tiffImg.open() % Open image
        
        % Iterate over y blocks
        for blck_y = 1:tiffImg.numBlcks(1)
            
            % Get image block
            [I,x,y] = getBlock(tiffImg,blck_x,blck_y);
            
            % Get threshold for block.
            if ~backgroundFirst || ~foregroundFirst
                threshold = tiffImg.threshold_fun(x,y);
            end

            if tiffImg.Use_GPU
                I = gpuArray(I);
            end
            
            % Smooth image
            Is = imfilter(I, tiffImg.Image_Smooth_Kernel, 'symmetric'); clear I
            
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
            
            % Take arrays off gpu before passing to partitioner
            if tiffImg.Use_GPU
                I = gather(Is); clear Is
                BWfull = gather(BW); clear BW
            else
                I = Is; clear Is
                BWfull = BW; clear BW
            end
            
            % Partition objects
            BW_part = partitioner.partition(BWfull, I, [tiffImg.xEdges(blck_x),tiffImg.yEdges(blck_y)]);
            writeBlock(tiffImg, maskID, logical(BW_part), blck_x, blck_y);            
        end % y block
        
        tiffImg.close(); % Prevent memory buildup
        progress.iteration_end(); % Update progress counter
    end % x block

    
catch ME
    tiffImg.close();
    tifflib('close',maskID)
    if tiffImg.Use_GPU
        gpuDevice([]);
    end
    fprintf('entered catch statement\n')
    rethrow(ME)
end

tiffImg.close();
tifflib('close',maskID)

end


function writeBlock(obj, maskID, I_blck, blck_x, blck_y)
tmpB_x_inds = obj.blck_x_inds + (blck_x-1)*obj.tilesPerBlck(2);
tmpB_y_inds = obj.blck_y_inds + (blck_y-1)*obj.tilesPerBlck(1);

if blck_x == obj.numBlcks(2)
    tmpB_x_inds(tmpB_x_inds > obj.numTiles(2)) = [];
end
if blck_y == obj.numBlcks(1)
    tmpB_y_inds(tmpB_y_inds > obj.numTiles(1)) = [];
end

cTiles = obj.tiles(tmpB_y_inds,tmpB_x_inds);
szcTiles = size(cTiles);

szI = size(I_blck);

for t = 1:numel(cTiles)
    
    [j,k] = ind2sub(szcTiles,t);
    
    tmpT_y_inds = obj.tile_y_inds + (j-1)*obj.tileSize(1);
    tmpT_x_inds = obj.tile_x_inds + (k-1)*obj.tileSize(2);
    
    if tmpT_x_inds(end) > szI(2)
        tmpT_x_inds(tmpT_x_inds > szI(2)) = [];
    end
    if tmpT_y_inds(end) > szI(1)
        tmpT_y_inds(tmpT_y_inds > szI(1)) = [];
    end
        
    tifflib('writeEncodedTile', maskID, cTiles(t), I_blck(tmpT_y_inds,tmpT_x_inds));
end
end