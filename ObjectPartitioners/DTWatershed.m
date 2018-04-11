classdef DTWatershed < ObjectPartitioner
    properties
        Smoothing_Size(1,1) double
    end
    properties (Hidden)
        disk
        H
    end
    methods
        function obj = DTWatershed(options)

            % Set defaults
            obj.Channel = "DAPI";
            obj.Smoothing_Size = 4;
%             obj.disk = strel('disk',2);
            obj.H = fspecial('gaussian',7*obj.Smoothing_Size,obj.Smoothing_Size);
            
            % Overwrite with user supplied options.
            if nargin > 0
                fd = fieldnames(options);

                for i = 1:numel(fd)
                    switch fd{i}
                        case 'Smoothing_Size'
                            obj.Smoothing_Size = options.Smoothing_Size;
                        case 'Area_Normalizer'
                            obj.Area_Normalizer = options.Area_Normalizer;
                            try
                                obj.Area_Normalizer(1,1);
                            catch
                                error('SALRGeo:Bad_area_normalizer', 'The Area_Normalizer should be a function that accepts two inputs: the x and y location of an object.')
                            end
                        case 'Intensity_Normalizer'
                            obj.Intensity_Normalizer = options.Intensity_Normalizer;
                            try
                                obj.Area_Normalizer(1,1);
                            catch
                                error('SALRGeo:Bad_intensity_normalizer', 'The Intensity_Normalizer should be a function that accepts two inputs: the x and y location of an object.')
                            end
                        case 'Restricted_Partitioning'
                            obj.Restricted_Partitioning = options.Restricted_Partitioning;
                        case 'Channel'
                            obj.Channel = options.Channel;
                    end
                end
            end
        end
        
        function obj = set.Smoothing_Size(obj, ss)
            obj.Smoothing_Size = ss;
            obj.disk = strel('disk',ss); %#ok<MCSUP>
%             obj.H = fspecial('gaussian',7*obj.Smoothing_Size,obj.Smoothing_Size); %#ok<MCSUP>
        end
        
        function BW_partitioned = partition(obj, BW, I, Image_Offset)
            % DTWATERSHED.PARTITION Partition overlapping nuclei by
            % locaing nuclei centers using the distance transform and
            % watershed.
            %
            % BW_partitioned = DTWatershed.partition(BW, I, Image_Offset)
            %
            % Input 
            %   BW : initial mask with partitially overlapping objects
            %   I : image
            %   Image_Offset : (optional) the offset of the image, if it
            %     comes from a larger image.
            %
            % Output
            %   BW_partitioned : The partitioned object mask
            
            % James Kapaldo
            
            if nargin < 3
                Image_Offset = [0 0];
            end
            
            % Pre-process the mask
            [CC, attempt] = preProcess(obj, I, BW, Image_Offset);

            % Pixel list of objects to attempt
            CC_attempt = CC.PixelIdxList(attempt);
            se = obj.disk;
%             h = obj.H;
            nRows = size(BW,1);
            
            if obj.Use_Parallel
                parfor i = 1:numel(CC_attempt)
                    BWi = createObjectImages(CC_attempt{i},nRows,true(numel(CC_attempt{i}),1));
                    D = bwdist(~BWi);
%                     D = -imfilter(D,h,'symmetric');
                    D = -imopen(D,se);
                    D(~BWi) = inf;
                    L = watershed(D,8);
                    BWi_p = (L>0) & BWi;
                    CC_attempt{i} = CC_attempt{i}(BWi_p(BWi)); 
                end
            else
                for i = 1:numel(CC_attempt)
                    BWi = createObjectImages(CC_attempt{i},nRows,true(numel(CC_attempt{i}),1));
                    D = bwdist(~BWi);
%                     D = -imfilter(D,h,'symmetric');
                    D = -imopen(D,se);
                    D(~BWi) = inf;
                    L = watershed(D,8);
                    BWi_p = (L>0) & BWi;
                    CC_attempt{i} = CC_attempt{i}(BWi_p(BWi)); 
                end
            end
            
            BW_partitioned = false(size(BW));
            BW_partitioned(cat(1,CC_attempt{:})) = true;
% 
%             % partition objects
%             D = bwdist(~BW_a);
%             D = -imopen(D,obj.disk);
%             D(~BW_a) = inf;
%             L = watershed(D,8);
%             BW_partitioned = (L > 0) & BW_a;
            
            % add objects that we did not need to attempt partitioning back
            % in.
            CC_good = CC.PixelIdxList(~attempt);
            BW_partitioned(cat(1,CC_good{:})) = true;

            % Post-process the mask
            BW_partitioned = postProcess(obj, BW_partitioned);
        end
    end
end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
