classdef (Abstract) FeatureGroup < matlab.mixin.Heterogeneous
    properties
        Channel(1,1) string = ""
        Options(1,1) struct
    end
    properties (SetAccess = protected)
        GroupName(1,1) string
    end
    properties (SetAccess = protected, Hidden)
        requiredOptions(1,:) cell
    end
    properties (Abstract, Dependent, SetAccess = protected)
        FeatureNames
    end
    
    methods (Abstract)
        x = Compute(obj);
    end
    
    methods
        function obj = set.Options(obj, options)
            requirementsCheck(obj, options)
            for i = 1:numel(obj.requiredOptions) %#ok<*MCSUP>
                obj.Options.(obj.requiredOptions{i}) = options.(obj.requiredOptions{i});
            end
        end
    end
    
    methods (Access=protected)
        function requirementsCheck(obj, options)
            
            fd = fieldnames(options)';
            
            for i = 1:numel(obj.requiredOptions)
                if ~contains(fd, obj.requiredOptions{i})
                    error('FeatureGroup:missingOptions','Missing options. %s FeatureGroup requires an option named, %s.', obj.GroupName, obj.requiredOptions{i})
                end
            end
            
%             if numel(fd) > numel(obj.requiredOptions)
%                 warning('FeatureGroup:extraOptions','Extra options given for FeatureGroup, %s. Extra options will be ignored.', obj.GroupName);
%             end
        end
        
%         function reload(obj,S)
%             obj.Channel = S
%         end
    end
    
    methods (Static, Sealed, Access = protected)
        function default_object = getDefaultScalarElement
            default_object = EmptyFeature;
        end
    end
    
    methods (Sealed)
        function tf = isempty(obj)
            tf = false(size(obj));
            for i = 1:numel(obj) 
                tf(i) = isa(obj(i),'EmptyFeature');
            end
        end
        
        function S = saveobj(obj)
%             warning('off','MATLAB:structOnObject');
%             S = struct(obj);
            S.Channel = obj.Channel;
            S.GroupName = obj.GroupName;
            S.Options = obj.Options;
            S.requiredOptions = obj.requiredOptions;
%             warning('on','MATLAB:structOnObject');
        end
    end
    
    methods (Static)
        function obj = loadobj(S)
            if isstruct(S)
                obj = eval(S.GroupName);
                obj.Channel = S.Channel;
                obj.Options = S.Options;
                obj.requiredOptions = S.requiredOptions;
            else
                obj = S;
            end
            
        end
    end
    
    
end
%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
