classdef SALRGeoPartition < ObjectPartitioner
    properties
        % Options for SALR particle clustering
        SALR_Options seedPointOptions

        % Options for geometric partitioning
        Part_Options partitionOptions
    end

    methods
        function obj = SALRGeoPartition(options)

            % Set default options.
            obj.Channel = "DAPI";
            % SALR clustering options
            try
                obj.SALR_Options = seedPointOptions();
            catch
                error('SALRGeo:missingDependance','SALRGeo requires the SALR clustering code.');
            end
            obj.SALR_Options.Wigner_Seitz_Radius = 5;
            obj.SALR_Options.Point_Selection_Method = 'r0set_uniformRandom';
            obj.SALR_Options.Maximum_Initial_Potential = 1/5;
            obj.SALR_Options.Max_Distance_Transform = 18;
            obj.SALR_Options.Potential_Parameters = [-1, 2, 13];
            obj.SALR_Options.Use_Parallel = true;

            % Geometric partitioning options
            try
                obj.Part_Options = partitionOptions();
            catch
                error('SALRGeo:missingDependance','SALRGeo requires the geometric partitioning code.');
            end
            obj.Part_Options.Use_GPU = 1;
            obj.Part_Options.Minimum_Hole_Size = 10;
            obj.Part_Options.Minimum_Object_Size = 150;
                
            % Overwrite with user supplied options.
            if nargin > 0
                fd = fieldnames(options);

                for i = 1:numel(fd)
                    switch fd{i}
                        case 'SALR_Options'
                            obj.SALR_Options = options.SALR_Options;
                        case 'Part_Options'
                            obj.Part_Options = options.Part_Options;
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
            
            obj.Minimum_Object_Size = obj.Part_Options.Minimum_Object_Size;
            obj.Minimum_Hole_Size = obj.Part_Options.Minimum_Hole_Size;
            obj.Use_Parallel = obj.SALR_Options.Use_Parallel;
        end

        function BW_partitioned = partition(obj, BW, I, Image_Offset)
            % SALRGEOPARTITION.PARTITION Partition overlapping nuclei by
            % locaing nuclei centers using SALR particle clustering and
            % then applying a seed-point based geometric partition
            %
            % BW_partitioned = SALRGeoPartition.partition(BW, I, Image_Offset)
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
            
            if nargin < 4
                Image_Offset = [0 0];
            end
            
            % Pre-process the mask
            [CC, attempt] = preProcess(obj, I, BW, Image_Offset);
            
            % Pixel list of objects to attempt
            CC_attempt = CC;
            CC_attempt.PixelIdxList = CC_attempt.PixelIdxList(attempt);
            CC_attempt.NumObjects = sum(attempt);
            
            
            if CC_attempt.NumObjects
                % Partition the objects
                BW_partitioned = declumpNuclei(I, CC_attempt, obj.SALR_Options, obj.Part_Options);%, attempt);
            else
                BW_partitioned = false(size(BW));
            end

            BW_partitioned(cat(1,CC.PixelIdxList{~attempt})) = true;
            
            % Post-process the mask
            BW_partitioned = postProcess(obj, BW_partitioned);
        end

        
    end
end

%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
