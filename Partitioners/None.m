classdef None < ObjectPartitioner
    methods (Static)
        function trueFalse = attemptPartitioning(objInfo)
            trueFalse = true(size(objInfo,1),1);
        end
    end
    methods
        function obj = None()
        end
        
        function mask_partitioned = partition(obj, mask, ~)
            mask_partitioned = mask;
        end
    end
end