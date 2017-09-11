classdef FeatureExtractor
    properties
        objectPartioner ObjectPartioner = None
        featureGroups(1,:) FeatureGroup
        channels(1,:) cell
    end
    
    methods
        function obj = FeatureExtractor(channels)
            obj.channels = channels;
        end
        
        function obj = addFeatureGroup(obj,group)
%             if ~isa(group, 'FeatureGroup')
%                 error('FeatureExtractor:badGroup','featureGroups must be a subclass of FeatureGroup');
%             end
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
            L = bwlabel(mask);
            
            % Extract features and feature names.
            x = [];
            names = [];
            
            for i = 1:numel(obj.featureGroups)
                ch = find(strcmp(obj.channels, obj.featureGroups(i).Channel));
                if isempty(ch)
                    xi = obj.featureGroups(i).Compute([], L);
                else
                    xi = obj.featureGroups(i).Compute(I(:,:,ch), L);
                end
                x = [x, xi];
                names = [names, obj.featureGroups(i).FeatureNames];
            end
        end
    end 
end