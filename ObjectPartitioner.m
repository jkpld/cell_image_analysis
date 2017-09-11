classdef (Abstract) ObjectPartitioner
    properties
        Options
        Channel = ''
    end
    
    methods (Abstract, Static)
        trueFalse = attemptPartioning(obj)
    end
    
    methods (Abstract)
        mask_partioned = partition(obj)
    end
end