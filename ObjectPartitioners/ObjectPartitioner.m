classdef (Abstract) ObjectPartitioner < matlab.mixin.Heterogeneous
    properties
        Channel(1,1) string = ""
        
        % Area_Normalizer : function that takes in image locations and
        % outputs the median object area at those positions. The areas are
        % normalized before going through the attemptPartitioning()
        % functions.
        Area_Normalizer = [];

        % Intensity_Normalizer : function that takes in image locations and
        % outputs the median object intensity at those positions. The
        % intensities are normalized before going through the
        % attemptPartitioning() functions.
        Intensity_Normalizer = [];
        
        % Restricted_Partitioning : Logical flag. If true, then the object
        % area and intensitiy will be normalized (using the Area_Normalizer
        % and Intensity_Normalizer) and the results passed to the
        % attemptPartitioning() method. If false, then all objects will be
        % kept, and all partitioning will be attempted on all objects.
        Restricted_Partitioning(1,1) logical = false;
        
        Minimum_Object_Size(1,1) double = 150;
        Minimum_Hole_Size(1,1) double = 10;
        Use_Parallel(1,1) logical = false;
        
        % Acceptance_Contour : Nx2 array descripting a contour in the
        % normalized DAPI-Area plane. Objects inside the contour will
        % attempt to be partitioned. Objects outside of the contour will be
        % ignored.
        Acceptance_Contour = [3, 1.1; 10, 4; 10 7; 5 7; 2 3; 1.6 2.5; 1.6 1.1; 3 1.1];
    end
    
    methods (Abstract)
        mask_partioned = partition(obj, BW, I, Image_Offset)
    end
    
    methods
        function tf = attemptPartitioning(obj, A, I)
            % ATTEMPTPARTITIONING Return true if the objects' normalized
            % area and intensity are inside the Acceptance_Contour
            %
            % This function is only valid when the nuclei data is
            % normalized such that the G1 nuclei have an area of 1 and an
            % intensity of 1.
            
            tf = inpolygon(I,A,obj.Acceptance_Contour(:,1),obj.Acceptance_Contour(:,2));
        end
        
        function [CC, attempt] = preProcess(obj, I, BW, Image_Offset)
            % PREPROCESS Determine the objects to try and partition and
            % remove small holes from the objects.
            %
            % [CC, attempt] = preProcess(obj, I, BW, Image_Offset)
            
            if nargin < 3
                Image_Offset = [0 0];
            end
            
            nRows = size(BW,1); % Image size
            
            % Remove small holes from the mask
            BW = ~bwareaopen(~BW, obj.Minimum_Hole_Size, 4);
            
            % Connected components
            CC = bwconncomp(BW);
            
            if obj.Restricted_Partitioning
                pixIdx = CC.PixelIdxList;

                % Object area, centroid, and intensity
                A = zeros(CC.NumObjects,1);
                mu = zeros(CC.NumObjects,2);
                iI = zeros(CC.NumObjects,1);

                if obj.Use_Parallel
                    Is = cellfun(@(x) I(x), pixIdx, 'UniformOutput',false);
                    parfor i = 1:numel(pixIdx)
                        A(i) = numel(pixIdx{i});
                        y = rem(pixIdx{i} - 1, nRows) + 1;
                        x = (pixIdx{i} - y)/nRows + 1;
                        mu(i,:) = mean([x,y],1);
                        iI(i) = sum(Is{i});
                    end
                else
                    for i = 1:numel(pixIdx)
                        A(i) = numel(pixIdx{i});
                        y = rem(pixIdx{i} - 1, nRows) + 1;
                        x = (pixIdx{i} - y)/nRows + 1;
                        mu(i,:) = mean([x,y],1);
                        iI(i) = sum(I(pixIdx{i}));
                    end
                end

                tooSmall = A < obj.Minimum_Object_Size;
                
                % Offset centroids by the image offset.
                mu = mu + Image_Offset;
                
                % Normalize area and intensity
                if ~isempty(obj.Area_Normalizer)
                    A = obj.Area_Normalizer(A, mu(:,1),mu(:,2));
                end
                if ~isempty(obj.Intensity_Normalizer)
                    iI = obj.Intensity_Normalizer(iI, mu(:,1),mu(:,2));
                end

                % Determine objects to attempt partitioning on
                attempt = attemptPartitioning(obj, A, iI);
            else                
                attempt = true(CC.NumObjects,1);
                tooSmall = cellfun(@numel,CC.PixelIdxList) < obj.Minimum_Object_Size;
            end
            
            attempt = attempt(~tooSmall);
            CC.PixelIdxList = CC.PixelIdxList(~tooSmall);
            CC.NumObjects = CC.NumObjects - sum(tooSmall);
        end
        
        function BW_part = postProcess(obj, BW)
            % Remove boundary objects and objects smaller than the minimum
            % allowed size.
            BW_part = bwareaopen(clear2dborder(BW), obj.Minimum_Object_Size);
        end
    end
    
    methods (Static, Sealed, Access = protected)
        function default_object = getDefaultScalarElement
            default_object = None;
        end
    end
    
end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
