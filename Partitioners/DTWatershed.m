classdef DTWatershed < ObjectPartitioner
    methods (Static)
        function trueFalse = attemptPartitioning(objInfo)
            trueFalse = true(size(objInfo,1),1);
        end
    end
    methods
        function obj = DTWatershed()
            obj.Options = struct(...
                'SmoothingSize',5, ... 
                'Minimum_Object_Size', 150);
        end
        
        function mask_partitioned = partition(obj, mask, ~)
            
            s = obj.Options.SmoothingSize;

            D = bwdist(~mask);
            D = -imopen(D,strel('disk',s));
            D(~mask) = inf;
            L = watershed(D,8);
            mask_partitioned = (L > 0) & BW;
                       
            mask_partitioned = bwareaopen(mask_partitioned, obj.Options.Minimum_Object_Size);
        end
    end
end