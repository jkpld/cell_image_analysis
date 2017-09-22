classdef None < ObjectPartitioner
    methods
        function obj = None()
        end
        
        function mask_partitioned = partition(obj, mask, ~)
            mask_partitioned = mask;
        end
    end
end
