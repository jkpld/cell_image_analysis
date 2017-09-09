classdef featureOptions
    properties
        features        = {};
        options         = struct();
    end
    properties (Hidden)
        calculateOutside = {};
        groupFeatureNames = {};
        channelNames = {};
    end
%     properties (Dependent, SetAccess = private, Hidden)
%         featureNames
%     end
    
    methods
        function obj = featureOptions(varargin)
            if nargin > 0
                if numel(varargin) == 1 && isstruct(varargin{1})
                    op = varargin{1};
                    fd = fieldnames(op)';
                    for fieldname = fd
                        obj = addFeature(obj, char(fieldname), op.(char(fieldname)));
                    end
                else
                    error('featureOptions:badInput','Input should be a structure of structures where the first level field names are the names of the features and the second level fields give the options for the features.')
                end
            end
        end

        function obj = addFeature(obj, featureName, options)
            switch featureName
                case 'Centroid'
                    if ~any(strcmp(obj.features, 'Centroid'))
                        obj.features = [obj.features, 'Centroid'];
                    end
                    obj.groupFeatureNames.Centroid = {'Centroid_x','Centroid_y'};
                case 'Area'
                    if ~any(strcmp(obj.features, 'Area'))
                        obj.features = [obj.features, 'Area'];
                    end
                    obj.groupFeatureNames.Area = {'Shape_Area'};
                case 'Shape'
                    if ~any(strcmp(obj.features, 'Shape'))
                        obj.features = [obj.features, 'Shape'];
                    end
                    if ~any(strcmp(obj.calculateOutside, 'Boundary'))
                        obj.calculateOutside = [obj.calculateOutside, 'Boundary'];
                    end
                    if ~any(strcmp(obj.calculateOutside, 'FourierDescriptors'))
                        obj.calculateOutside = [obj.calculateOutside, 'FourierDescriptors'];
                    end
                    
                    fd = fieldnames(options)';
                    requiredFields = {'Curvature_Smoothing_Size', 'Number_Boundary_Vertices','Number_Fourier_Descriptors'};
                    for i = 1:numel(requiredFields)
                        if ~any(strcmp(fd,requiredFields{i}))
                            error('featureOptions:missingFields','Mission options. Shape feature requires an option named, %s.', requiredFields{i})
                        end
                    end
                    
                    % Create filters for getting derivatives
                    k = options.Curvature_Smoothing_Size;

                    filtSize = round(7*k);
                    l = -floor(filtSize/2):ceil(filtSize/2);
                    G = exp(-(l).^2/(2*k^2))/(k*sqrt(2*pi));
                    dG = -(l) .* G / k^2;
                    ddG = - G / k^2 + (l).^2 .* G / k^4;

                    G = G(:);
                    dG = dG(:);
                    ddG = ddG(:);
                    
                    obj.options.G = G;
                    obj.options.dG = dG;
                    obj.options.ddG = ddG;
                    obj.options.FD_Number_Boundary_Vertices = options.Number_Boundary_Vertices;
                    obj.options.FD_Number_Fourier_Descriptors = options.Number_Fourier_Descriptors;
                    
                    obj.groupFeatureNames.Shape = cellstr('Shape_' + [string({'MajorAxisLength','MinorAxisLength','Eccentricity','Orientation','MeanPositiveCurvature','MeanNegativeCurvature','Solidity','FormFactor'}), 'FourierDescriptor_'+string(1:options.Number_Fourier_Descriptors)]);
                    
                case 'Intensity'
                    if ~any(strcmp(obj.features, 'Intensity'))
                        obj.features = [obj.features, 'Intensity'];
                    end
                    if ~any(strcmp(obj.calculateOutside, 'SlicedIntensity'))
                        obj.calculateOutside = [obj.calculateOutside, 'SlicedIntensity'];
                    end
                    
                    obj.groupFeatureNames.Intensity = cellstr('Intensity_' + string({'Min', 'Mean', 'Max', '1stQuartile','Median','3rdQuartile','Std','MAD'}));
                    
                case 'RadialIntensity'
                case 'Texture'
                case 'Granularity'
                case 'LBP'
                case 'Correlation'
                case 'SimplifiedBoundary'
                    if ~any(strcmp(obj.features, 'SimplifiedBoundary'))
                        obj.features = [obj.features, 'SimplifiedBoundary'];
                    end
                    if ~any(strcmp(obj.calculateOutside, 'Boundary'))
                        obj.calculateOutside = [obj.calculateOutside, 'Boundary'];
                    end
                    if ~any(strcmp(obj.calculateOutside, 'SimplifiedBoundary'))
                        obj.calculateOutside = [obj.calculateOutside, 'SimplifiedBoundary'];
                    end
                    
                    fd = fieldnames(options)';
                    requiredFields = {'Curvature_Smoothing_Size', 'Number_Boundary_Vertices','Number_Fourier_Descriptors'};
                    for i = 1:numel(requiredFields)
                        if ~any(strcmp(fd,requiredFields{i}))
                            error('featureOptions:missingFields','Mission options. SimplifiedBoundary feature requires an option named, %s.', requiredFields{i})
                        end
                    end
                    
                    obj.options.SB_Number_Boundary_Vertices = options.Number_Boundary_Vertices;
                    obj.options.SB_Number_Fourier_Descriptors = options.Number_Fourier_Descriptors;
                    
                    obj.groupFeaturenames.SimplifiedBoundary = cellstr('SimplifiedBoundary_FourierDescriptor_'+string(1:options.Number_Fourier_Descriptors));
                    
                case 'BoundingBox'
                    if ~any(strcmp(obj.features, 'BoundingBox'))
                        obj.features = [obj.features, 'BoundingBox'];
                    end
                    obj.groupFeaturenames.BoundingBox = cellstr('BoundingBox_' + string({'x','y','w','h'}));
                otherwise
                    warning('featureOptions:unknownFeature','Unknown feature, %s. Feature is being ignored.', featureName)
            end
            
        end % addFeature
        
        function obj = set.channelNames(obj,channelNames)
            
            obj.channelNames = channelNames;
                
        end % setChannelNames
        
        function featureNameList = featureNames(obj)
            
            featureNameList = {};
            for i = 1:numel(obj.features)
                switch obj.features{i}
                    case 'Centroid'
                        featureNameList = [featureNameList, obj.groupFeatureNames.Centroid]; %#ok<*AGROW>
                    case 'Area'
                        featureNameList = [featureNameList, obj.groupFeatureNames.Area];
                    case 'Shape'
                        featureNameList = [featureNameList, obj.groupFeatureNames.Shape];
                    case 'Intensity'
                        if ~isempty(obj.channelNames)
                            intensityFeatureNames = (string(obj.groupFeatureNames.Intensity) + '_' + string(obj.channelNames)')';
                            featureNameList = [featureNameList, cellstr(intensityFeatureNames(:)')];
                        else
                            featureNameList = [featureNameList, obj.groupFeatureNames.Intensity];
                        end
                    case 'RadialIntensity'
                    case 'Texture'
                    case 'Granularity'
                    case 'LBP'
                    case 'Correlation'
                    case 'SimplifiedBoundary'
                        featureNameList = [featureNameList, obj.groupFeatureNames.SimplifiedBoundary];
                    case 'BoundingBox'
                        featureNameList = [featureNameList, obj.groupFeatureNames.BoundingBox];
                end
            end
            
        end
    end
    
end