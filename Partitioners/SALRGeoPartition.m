classdef SALRGeoPartition < ObjectPartitioner
    properties
        % Options for SALR particle clustering
        SALR_Options seedPointOptions

        % Options for geometric partitioning
        Part_Options partitionOptions

        % Area_Normalizer : function that takes in image locations and
        % outputs the median object area at those positions. The areas are
        % normalized before going through the keepObject() and
        % attemptPartitioning() functions.
        Area_Normalizer

        % Intensity_Normalizer : function that takes in image locations and
        % outputs the median object intensity at those positions. The
        % intensities are normalized before going through the keepObject()
        % and attemptPartitioning() functions.
        Intensity_Normalizer
        
        % Restricted_Partitioning : Logical flag. If true, then the object
        % area and intensitiy will be normalized (using the Area_Normalizer
        % and Intensity_Normalizer) and the results passed to the methods
        % keepObject() and attemptPartitioning(). If false, then all
        % objects will be kept, and all partitioning will be attempted on
        % all objects.
        Restricted_Partitioning(1,1) logical
    end

    methods
        function obj = SALRGeoPartition(options)

            % Set default options.
            
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

            obj.Area_Normalizer = @(~,~) 1;
            obj.Intensity_Normalizer = @(~,~) 1;
            
            obj.Restricted_Partitioning = false;
                
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
                    end
                end
            end
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
            
            if nargin < 3
                Image_Offset = [0 0];
            end
            
            nRows = size(BW,1); % Image size
            
            % Remove small holes from the mask
            BW = ~bwareaopen(~BW, obj.Part_Options.Minimum_Hole_Size, 4);
            
            % Connected components
            CC = bwconncomp(BW);
            
            if obj.Restricted_Partitioning
                Use_Parallel = obj.SALR_Options.Use_Parallel;
                pixIdx = CC.PixelIdxList;

                % Object area, centroid, and intensity
                A = zeros(CC.NumObjects,1);
                mu = zeros(CC.NumObjects,1);
                iI = zeros(CC.NumObjects,1);

                if Use_Parallel
                    Is = cellfun(@(x) I(x), pixIdx, 'UniformOutput',false);
                    parfor i = 1:numel(pixIdx)
                        A(i) = numel(pixIdx{i});
                        y = rem(pixIdx{i} - 1, nRows) + 1;
                        x = (pixIdx{i} - y)/nRows + 1;
                        mu(i,:) = mean([x,y]);
                        iI(i) = sum(Is{i});
                    end
                else
                    for i = 1:numel(pixIdx)
                        A(i) = numel(pixIdx{i});
                        y = rem(pixIdx{i} - 1, nRows) + 1;
                        x = (pixIdx{i} - y)/nRows + 1;
                        mu(i,:) = mean([x,y]);
                        iI(i) = sum(I(pixIdx{i}));
                    end
                end

                % Offset centroids by the image offset.
                mu = mu + Image_Offset;

                % Objects that are too small
                keep = A > obj.Part_Options.Minimum_Object_Size;
                
                % Normalize area and intensity
                A = A ./ obj.Area_Normalizer(mu(:,1),mu(:,2));
                iI = iI ./ obj.Intensity_Normalizer(mu(:,1),mu(:,2));

                % Determine objects that should be kept
                keep = keep & keepObject(obj, A, iI);
                
                A = A(keep);
                iI = iI(keep);
                CC.PixelIdxList = CC.PixelIdxList(keep);
                CC.NumObjects = sum(keep);

                % Determine objects to attempt partitioning on
                attempt = attemptPartitioning(obj, A, iI);
            else
                A = cellfun(@numel, CC.PixelIdxList);
                keep = A > obj.Part_Options.Minimum_Object_Size;
                CC.PixelIdxList = CC.PixelIdxList(keep);
                CC.NumObjects = sum(keep);
                
                attempt = true(CC.NumObjects,1);
            end

            % Partition the objects
            BW_partitioned = declumpNuclei(I, CC, obj.SALR_Options, obj.Part_Options, attempt);

            % Remove boundary objects and objects smaller than the minimum
            % allowed size.
            BW_partitioned = bwareaopen(imclearborder(BW_partitioned), obj.Part_Options.Minimum_Object_Size);
        end

        function tf = attemptPartitioning(~, A, I)
            tf = (A > 1.1) & (I > 1.6);
        end

        function tf = keepObject(~, A, I)
            acceptanceContour = [3, 1; 10, 5; 10 7; 8 7; 2 3; 1 2; 0 2; 0 0; 3 0; 3 1];
            tf =  inpolygon(I,A,acceptanceContour(:,1),acceptanceContour(:,2));
        end
    end
end
