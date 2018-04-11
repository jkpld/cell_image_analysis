classdef None < ObjectPartitioner
    methods
        function obj = None()
        end
        
        function mask_partitioned = partition(obj, mask, ~, ~)
            mask_partitioned = mask;
        end
    end
end

%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
