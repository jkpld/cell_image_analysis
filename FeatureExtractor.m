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
        
        function [x, names, mask] = Compute(obj, mask, I, Image_Offset)
            if size(I,3) ~= numel(obj.channels)
                error('FeatureExtractor:channelsDimensionMismatch','The input image I must have the same number of channels (and be in the same order) as the channels cell array used when creating the object.');
            end
            
            if nargin < 4
                Image_Offset = zeros(1,2,'single');
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
            
            if numel(groups) == 0
                x = [];
                names = [];
                warning('Compute:noValidGroups','There are no valid groups - either the FeatureGroups are EmptyFeatures or the FeatureGroup channel names to not match the FeatureExtractor channel names.')
                return;
            end
            
            % Get feature names
            names = {groups.FeatureNames};
            num_xi = cellfun(@numel, names);
            names = cat(2,names{:});
            
            % Initialize feature array
            x = zeros(max(L(:)), sum(num_xi), 'single');

            if size(x,1) == 0
                return;
            end
            
            % feature end indices
            num_xi = cumsum([0,num_xi]);
            
            for i = 1:numel(groups)
                ch = obj.channels == groups(i).Channel;
                if ~any(ch)
                    xi = groups(i).Compute([], L);
                else
                    xi = groups(i).Compute(I(:,:,ch), L);
                end

                if isa(groups(i),'Location')
                    xi(:,1:4) = xi(:,1:4) + [Image_Offset, Image_Offset];
                end
                if isa(groups(i),'BasicProps')
                    xi(:,1:2) = xi(:,1:2) + Image_Offset;
                end
                x(:, num_xi(i)+1:num_xi(i+1)) = xi;
            end
        end
    end 
end