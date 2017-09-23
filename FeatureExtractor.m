classdef FeatureExtractor
    properties
        objectPartioner ObjectPartitioner = None
        featureGroups(1,:) FeatureGroup
        channels(1,:) string
    end
    
    methods
        function obj = FeatureExtractor(channels)
            if nargin > 0
                obj.channels = channels;
            end
        end
        
        function obj = addFeatureGroup(obj,group)
            obj.featureGroups(end+1) = group;
        end
        
        function [x, names, mask] = Compute(obj, mask, I)
            if size(I,3) ~= numel(obj.channels)
                error('FeatureExtractor:channelsDimensionMismatch','The input image I must have the same number of channels (and be in the same order) as the channels cell array used when creating the object.');
            end
            
            % Partition the objects?
            if ~isa(obj.objectPartioner,'None')               
                ch = find(strcmp(obj.channels, obj.objectPartioner.Channel));
                if isempty(ch)
                    mask = obj.objectPartioner.partition(mask, []);
                else
                    mask = obj.objectPartioner.partition(mask, I(:,:,ch));
                end
            end
            
            % Create label matrix
            L = single(bwlabel(mask));
            
            if isinteger(I)
                I = single(I)/single(intmax(class(I)));
            else
                I = single(I);
            end
            
            % Remove empty groups
            validChannels = ["", obj.channels];
            groups = obj.featureGroups(~isempty(obj.featureGroups) & ismember([obj.featureGroups.Channel],validChannels));
            
            % Get feature names
            names = {groups.FeatureNames};
            num_xi = cellfun(@numel, names);
            
            % Initialize feature array
            x = zeros(max(L(:)), sum(num_xi), 'single');

            % feature end indices
            num_xi = cumsum([0,num_xi]);

            
            for i = 1:numel(groups)
                ch = obj.channels == groups(i).Channel;
                if ~any(ch)
                    xi = groups(i).Compute([], L);
                else
                    xi = groups(i).Compute(I(:,:,ch), L);
                end

                x(:, num_xi(i)+1:num_xi(i+1)) = xi;
            end
            
            names = cat(2,names{:});
        end
    end 
end