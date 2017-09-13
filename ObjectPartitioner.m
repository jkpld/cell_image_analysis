classdef (Abstract) ObjectPartitioner
    properties
        Channel(1,1) string = ""
    end
    
    methods (Abstract)
        mask_partioned = partition(obj, BW, I, Image_Offset)
        
        tf = keepObject(obj) % Take in information about each object (before partitioning) and determine if the object should be kept (true) or thrown away (false)
        tf = attemptPartitioning(obj) % Take in information about each object (before partitioning) and determine if partitioning should be attempted (E.g. if the object is smaller than some value, then you might know you do not need to partition in.)
    end
end