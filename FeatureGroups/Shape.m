classdef Shape < FeatureGroup
    % SHAPE Compute geometric features describing the shape of each object.
    % The features are
    %   Area - the object area
    %   Major axis length - the major axis length of an ellipse fit to the
    %     object.
    %   Minor axis length - the minor axis length of an ellipse fit to the
    %     object.
    %   Eccentricity - the eccentricity of the ellipse fit to the object.
    %   Solidity - the area of the object divided by the area of the
    %     object's convex hull
    %   Form factor - the area of the object divided by the area of a
    %     circle with the same circumfrence as the object's perimeter.
    %   Fourier descriptors - translation, rotation, and scale invarient
    %     Fourier descriptors that describe the object's shape.
    %
    % Properties :
    %
    % Channel - Not used. This FeatureGroup works directly on a mask.
    % requiredOptions - Structure with the following fields
    %   Number_Boundary_Vertices : The number of boundary vertices to use.
    %     All boundaries will be resized to this size before computing the
    %     Fourier descriptors and perimeters, which are used in the form
    %     factor feature.
    %   Number_Fourier_Descriptors : The number of Fourier descriptors to
    %     keep. (This must be less than Number_Boundary_Vertices.)
    %   Use_GPU : logical flag determining if a GPU is used to speed up the
    %     computation.
    %   Use_Parallel : logical flag determining if multiple CPUs will be
    %     used to speed up the computation. (This is used in computing the
    %     Solidity feature, since convex hull computation is not available
    %     on a GPU in Matlab.)
    %
    % Methods :
    %
    % Compute(~, L) - Compute features describing the shape of each object
    %   in the label matrix, L.
    %   x = Shape.Compute([], L)
    
    properties (Dependent, SetAccess = protected)
        FeatureNames
    end
    properties (Access = private, Hidden)
        Locator
    end
    methods
        function obj = Shape(channel,options)
            obj.GroupName = class(obj);
            
            % Default options
            defaultOptions = struct( ...
                'Number_Boundary_Vertices', 100, ...
                'Number_Fourier_Descriptors', 10, ...
                'Use_GPU', 0, ...
                'Use_Parallel', 0);
            obj.requiredOptions = fieldnames(defaultOptions)';
            
            if nargin < 2
                options = defaultOptions;
            end
            obj.Options = options;
            
            % Check the options to make sure they are valid.
            if obj.Options.Number_Fourier_Descriptors < 3
                warning('Shape:NFD_lt_3', 'The minimum number of Fourier descriptors is 3, because two of the Fourier descriptors are constrainted (and removed) to achive translation and scale invariance. The number of Fourier descriptors has been icnreased to 3.')
                obj.Options.Number_Fourier_Descriptors = 3;
            end
            if obj.Options.Number_Fourier_Descriptors < obj.Options.Number_Boundary_Vertices
                warning('Shape:NFD_lt_NBV', 'The number of Fourier descriptors must be less than the number of boundary vertices. The number of Fourier descriptors has been reduced.')
                obj.Options.Number_Fourier_Descriptors = obj.Options.Number_Boundary_Vertices;
            end
            
            if nargin >= 1
                if (channel ~= '') || ~isempty(channel)
                    warning('Shape:channelNotUsed', 'The Shape FeatureGroup does not operate on an image channel. Channel given will be ignored.')
                end
            end
            
            % Add a Location FeatureGroup that will be used to extract the
            % Location features if they are not provided to Compute.
            obj.Locator = Location([], struct('Use_GPU', obj.Options.Use_GPU));
            
        end
        
        function names = get.FeatureNames(obj)
            
            NFD = obj.Options.Number_Fourier_Descriptors;
            
            featList = [ ...
                "Area", ...
                "MajorAxisLength", ...
                "MinorAxisLength", ...
                "Eccentricity", ...
                "Solidity", ...
                "FormFactor", ...
                "FourierDescriptor_" + string(2:NFD-1)];
            
            names = obj.GroupName + "_" + featList;
        end
        
        function x = Compute(obj, ~, L)
            % SHAPE.COMPUTE Return features describing the shape of each
            % object
            %
            % x = Location.Compute([],L)
            % x = Location.Compute([],L,LocFeatures)
            %
            % Input
            %   L : label matrix
            %
            % Output
            %   x : N x (6 + Number_Fourier_Descriptors-2) array. N is the
            %     number of objects. The columns correspond to the features
            %     given by FeatureNames.
            %
            % Note : The number of Fourier descriptors returned is
            % Number_Fourier_Descriptors - 2, because two of the
            % Descriptors are used to make the features translation and
            % scale invariant, and provide no information.
            
            % James Kapaldo
                        
            % GPU, Parallel options
            useGPU = obj.Options.Use_GPU;
            useParallel = obj.Options.Use_Parallel;
            
            % ===========================================================
            % Object Boundaries
            %
            B = bwboundariesmex(double(L), 8);
            
            % ===========================================================
            % Fourier Descriptors
            %
            NBV = obj.Options.Number_Boundary_Vertices;
            NFD = obj.Options.Number_Fourier_Descriptors;
            
            FD = fourierDescriptors(B, NBV, NFD, useParallel);
            
            % ===========================================================
            % Pixel list & centroids
            %
            %   Note, this is a bit redundent if the Location FeatureGroup
            %   is used, but we need the object indices, and it doesn't
            %   take that long.
            
            % Initialize linear indices
            linIdx = single(1:numel(L))';
            
            % Get number of rows in image
            nRow = size(L,1);
            
            % Remove background pixels
            FG = L>0;
            L = single(L(FG));
            linIdx = linIdx(FG);
            
            % If using GPU, then send the arrays there
            if useGPU
                L = gpuArray(L);
                linIdx = gpuArray(linIdx);
                nRow = gpuArray(nRow);
            end
            
            N_obj = max(L); % number of objects
            
            % Row and column indices
            y = rem(linIdx-1, nRow) + 1;
            x = (linIdx - y)/nRow + 1; clear linIdx nRow;
                       
            % Pixels per object (also the area)
            if useGPU
                N = accumarray(L, gpuArray.ones(size(L),'single'), [N_obj, 1], @sum, single(1));
            else
                % Using sparse on CPU is faster than accumarray
                N = single(full(sparse(double(L),ones(size(L)), ones(size(L)), double(N_obj), 1)));
                N(N==0) = 1;
            end
            
            % Centroids
            x_mu = accumarray(L, x) ./ N;
            y_mu = accumarray(L, y) ./ N;

            
            % ===========================================================
            % Ellipse Parameters
            %
            %   Re-writen version from matlab's regionprops function.
            x = x - x_mu(L); clear x_mu
            y = y - y_mu(L); clear y_mu
            
            uxx = accumarray(L, x.^2) ./ N + 1/12;
            uyy = accumarray(L, y.^2) ./ N + 1/12;
            uxy = accumarray(L, x.*y) ./ N; clear x y
            
            tmp1 = sqrt( (uxx - uyy).^2 + 4*uxy.^2 ); clear uxy
            tmp2 = uxx + uyy; clear uxx uyy
            
            major = 2*sqrt(2)*sqrt(tmp2 + tmp1);
            minor = 2*sqrt(2)*sqrt(tmp2 - tmp1); clear tmp1 tmp2
            ecc = 2*sqrt( (major/2).^2 - (minor/2).^2 ) ./ major;
            
            if useGPU
                A = gather(N);
                e_param = [gather(major), gather(minor), gather(ecc)];
                clear N major minor eccentricity
            else
                A = N;
                e_param = [major, minor, ecc];
            end
            
            
            % ===========================================================
            % Solidity & Form factor
            %
            
            if useGPU
                N_obj = gather(N_obj); % Note, N_obj is still on GPU but cannot be accessed.
            end
            
            S_F = zeros(N_obj,2);
            if useParallel
                parfor i = 1:N_obj
                    S = ConvexArea(B{i});
                    F = computePerimeterFromBoundary(B{i});
                    S_F(i,:) = [S, F];
                end
            else
                for i = 1:N_obj
                    S = ConvexArea(B{i});
                    F = computePerimeterFromBoundary(B{i});
                    S_F(i,:) = [S, F];
                end
            end
            
            S_F(:,2) = S_F(:,2).^2 / (4*pi) ;
            S_F = A ./ S_F;
            
            % ===========================================================
            % Create final feature array.
            %
            x = [A, e_param, S_F, FD];
        end
    end
end

