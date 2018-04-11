classdef EmptyFeature < FeatureGroup
    % EMPTYFEATURE Empty feature group to act as a place holder in
    % heterogenous arrays of FeatureGroup's
    %
    % Properties :
    %
    % Channel - Not used.
    % requiredOptions - Empty structure
    %
    % Methods :
    %
    % Compute(~, ~) - Return an empty matrix
    %   x = EmptyFeature.Compute([], [])
    
    properties (Dependent, SetAccess = protected)
        FeatureNames
    end
    methods
        function obj = EmptyFeature(channel,options)
            obj.GroupName = class(obj);

            % Default options
            defaultOptions = struct();
            obj.requiredOptions = fieldnames(defaultOptions)';
            
            if nargin < 2
                options = defaultOptions;
            end
            obj.Options = options;
            
            if nargin >= 1
                if ~strcmp(channel, '') || ~isempty(channel)
                    warning('Location:channelNotUsed', 'The Location FeatureGroup does not operate on an image channel. Channel given will be ignored.')
                end
            end
        end
        
        function names = get.FeatureNames(~)           
            names = "";
        end
        
        function x = Compute(~, ~, ~)
            x = [];
        end
    end
end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
