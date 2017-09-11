classdef SALRGeo < ObjectPartitioner
    methods (Static)
        function trueFalse = attemptPartitioning(objInfo)
            trueFalse = true(size(objInfo,1),1);
        end
    end
    methods
        function obj = SALRGeo()
            % Geometric partitioning options
            partOptions = partitionOptions();
            partOptions.Use_GPU = 1;
            partOptions.Minimum_Hole_Size = 10;

            % SALR clustering options
            salrOptions = seedPointOptions();
            salrOptions.Wigner_Seitz_Radius = 5;
            salrOptions.Point_Selection_Method = 'r0set_uniformRandom';
            salrOptions.Maximum_Initial_Potential = 1/5;
            salrOptions.Max_Distance_Transform = 18;
            salrOptions.Potential_Parameters = [-1, 2, 13];
            salrOptions.Verbose = true;
            salrOptions.Use_Parallel = true;
            
            obj.Options.salrOptions = salrOptions;
            obj.Options.partOptions = partOptions;
            obj.Options.Minimum_Object_Size = 150;
            obj.Channel = 'DAPI';
        end
        
        function mask_partitioned = partition(obj, mask, I)
            
            Use_GPU = obj.Options.partOptions.Use_GPU;
            getBasicInfo = BasicProps('',struct('Use_GPU',Use_GPU));
            basicInfo = getBasicInfo.Compute(I, bwlabel(mask));
            
            attempt = attemptPartitioning(basicInfo);
            
            mask_partitioned = declumpNuclei(...
                I, ...
                mask, ...
                obj.Options.salrOptions, ...
                obj.Options.geoOptions, ...
                'tryPartitioning',attempt);
            
            mask_partitioned = bwareaopen(mask_partitioned, obj.Options.Minimum_Object_Size);
        end
    end
end