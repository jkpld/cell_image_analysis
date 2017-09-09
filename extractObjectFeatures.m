function features = extractObjectFeatures(CC,I,featureOptions, offset)



N = CC.NumObjects;
pixIdx = CC.PixelIdxList;

[m,n,ch] = size(I);

calculateOutside = featureOptions.calculateOutside;
features = featureOptions.features;
options = featureOptions.options;

SimplifiedBoundaryCalculated = false;
FourierDescriptorsCalculated = false;

boundariesCalculated = false;



for i = 1:numel(calculateOutside)
    switch calculateOutside{i}
        case 'Boundary'
            if ~boundariesCalculated
                L = labelmatrix(CC);
                B = bwboundariesmex(double(L), 8); % regionboundariesmex(double(L),8); 
                boundariesCalculated = true;
            end
            % I dont know what the difference between these two functions
            % are, but I want my boundaries to be clockwise oriented and I
            % used bwboundaries before with the fourier descriptors
            % function without having to check the orientation; so, i will
            % use that mex function.
        case 'FourierDescriptors'
            if ~boundariesCalculated
                L = labelmatrix(CC);
                B = bwboundariesmex(double(L), 8); % regionboundariesmex(double(L),8); 
                boundariesCalculated = true;
            end
            
            if ~FourierDescriptorsCalculated
                if any(strcmp(calculateOutside,'SimplifiedBoundary'))
                    if options.FD_Number_Boundary_Vertices == options.SB_Number_Boundary_Vertices && ...
                            options.FD_Number_Fourier_Descriptors == options.SB_Number_Fourier_Descriptors
                        SimplifiedBoundaryCalculated = true;
                        [FD, slicedSimplifiedBoundary] = fourierDescriptors(B,options.FD_Number_Boundary_Vertices,options.FD_Number_Fourier_Descriptors); %#ok<*ASGLU>
                    else
                        FD = fourierDescriptors(B,options.FD_Number_Boundary_Vertices,options.FD_Number_Fourier_Descriptors);
                    end
                else
                    FD = fourierDescriptors(B,options.FD_Number_Boundary_Vertices,options.FD_Number_Fourier_Descriptors);
                end
                FourierDescriptorsCalculated = true;
            end
        case 'SlicedIntensity'
            % Get the intensities of each channel for each object
            s = [];
            for k = 1:ch
                s = [s,sprintf('I(x+%d*n*m),',k-1)]; %#ok<AGROW>
            end
            s(end) = [];
            s = ['I = cellfun(@(x) cat(3,' s, '),pixIdx,''UniformOutput'',0);']; %#ok<AGROW>
            eval(s);
        case 'SimplifiedBoundary'
            if ~boundariesCalculated
                L = labelmatrix(CC);
                B = bwboundariesmex(double(L), 8); % regionboundariesmex(double(L),8); 
                boundariesCalculated = true;
            end
            
            if ~SimplifiedBoundaryCalculated
                if any(strcmp(calculateOutside,'FourierDescriptors'))
                    if options.FD_Number_Boundary_Vertices == options.SB_Number_Boundary_Vertices && ...
                            options.FD_Number_Fourier_Descriptors == options.SB_Number_Fourier_Descriptors
                        FourierDescriptorsCalculated = true;
                        [FD, slicedSimplifiedBoundary] = fourierDescriptors(B,options.FD_Number_Boundary_Vertices,options.FD_Number_Fourier_Descriptors);
                    else
                        [~, slicedSimplifiedBoundary] = fourierDescriptors(B,options.SB_Number_Boundary_Vertices,options.SB_Number_Fourier_Descriptors);
                    end
                else
                    [~, slicedSimplifiedBoundary] = fourierDescriptors(B,options.SB_Number_Boundary_Vertices,options.SB_Number_Fourier_Descriptors);
                end
                SimplifiedBoundaryCalculated = true;
            end
        otherwise
    end
end

if ~any(strcmp(calculateOutside,'Boundary'))
    B = cell(N,1);
end

% Initialize sliced variables
for i = 1:numel(features)
    switch features{i}
        case 'Centroid'
            slicedCentroid = nan(N,2);
        case 'Area'
            slicedArea = nan(N,1);
        case 'Shape'
            slicedShape = nan(N,8);
            
        case 'Intensity'
            slicedIntensity = nan(N,ch*8);
        case 'RadialIntensity'
        case 'Texture'
        case 'Granularity'
        case 'LBP'
        case 'Correlation'
        case 'BoundingBox'
            slicedBoundingBox = nan(N,4);
        otherwise
    end
end


% Start feature extraction

emptyIndependentFeatures = struct('Centroid',[], ...
                                  'Area', [], ...
                                  'BoundingBox',[], ...
                                  'PixelList',[], ...
                                  'LinearPixelList',[], ...
                                  'NumImageRows',m);
                      

parfor obj = 1:N
    
    objFeatures = emptyIndependentFeatures;
    objFeatures.LinearPixelList = pixIdx{obj};
    objFeatures.Area = numel(pixIdx{obj});

    for i = 1:numel(features)

        switch features{i}
            case 'Centroid'
                objFeatures = computeCentroid(objFeatures);                
                slicedCentroid(obj,:) = objFeatures.Centroid + offset; %#ok<*PFOUS>
            case 'Area' % Area
                
                slicedArea(obj) = objFeatures.Area;
                
            case 'Shape' % Ellipse parameters, Mean positive and negative curvature, Solidity, Form factor
                
                % Prerequisits: centroid, boundingbox
                objFeatures = computeCentroid(objFeatures);
                objFeatures = computeBoundingBox(objFeatures);
                
                % Compute ellipse parameters
                [mj,mn,ecc,ori] = computeEllipseParameters(objFeatures.PixelList,objFeatures.Centroid);
                
                % Compute boundary curvature means
                [meanKpos, meanKneg] = computeBoundaryCurvatureMeans(B{obj}, options);
                
                % Compute convex hull area
                K = convhull(B{obj}(:,1),B{obj}(:,2));
                hull = B{obj}(K,:);
                
                r = hull(:,2) - objFeatures.BoundingBox(2) + 1;
                c = hull(:,1) - objFeatures.BoundingBox(1) + 1;
                M = objFeatures.BoundingBox(4);
                N = objFeatures.BoundingBox(3);

                convexImage = roipoly(M,N,c,r);
                CHA = sum(convexImage(:));
                
                % Compute perimeter
                P = computePerimeterFromBoundary(B{obj});
                
                % Save the features
                slicedShape(obj,:) = [mj,mn,ecc,ori,meanKpos,meanKneg,objFeatures.Area/CHA,4*pi*objFeatures.Area/P^2]
                
            case 'Intensity' % min, mean, max, 1st quartile, median, 3rd quartile, std, mad
                minI = min(I{obj},[],1); % 1x1x3
                meanI = mean(I{obj},1); % 1x1x3
                maxI = max(I{obj},[],1); % 1x1x3
                prctileI = permute(prctile(I{obj},[25,50,75]),[2,1,3]); % 1x3x3
                stdI = std(I{obj},[],1); % 1x1x3
                madI = mad(I{obj},1); %median(I{obj}-prctileI(1,2,:),1); %  1x1x3

                objI = [minI, meanI, maxI, prctileI, stdI, madI]; % should be 1 x L x ch
                objI = permute(objI,[2,3,1]);

                slicedIntensity(obj,:) = objI(:)'; 
                
            case 'BoundingBox'
                objFeatures = computeBoundingBox(objFeatures);                
                slicedBoundingBox(obj,:) = objFeatures.BoundingBox + [offset 0 0];
        end
    end % for features    
end % parfor

if any(strcmp(features,'Shape'))
    slicedShape = [slicedShape, FD]; %#ok<NASGU>
end

% Concoctinate the slicedFeatures into one array
str = sprintf('sliced%s,',features{:});
str(end) = [];
eval(['features = [' str '];']);

end % function

function objFeatures = computeCentroid(objFeatures)
    if isempty(objFeatures.Centroid)
        objFeatures = computePixelList(objFeatures);
        objFeatures.Centroid = mean(objFeatures.PixelList);
    end
end

function objFeatures = computePixelList(objFeatures)
    if isempty(objFeatures.PixelList)
        pixelList(:,2) = rem(objFeatures.LinearPixelList-1, objFeatures.NumImageRows) + 1;
        pixelList(:,1) = (objFeatures.LinearPixelList - pixelList(:,2))/objFeatures.NumImageRows + 1;
        objFeatures.PixelList = pixelList;
    end
end

function objFeatures = computeBoundingBox(objFeatures)
    if isempty(objFatures.BoundingBox)
        objFeatures = computePixelList(objFeatures);

        topLeft = min(objFeatures.PixelList);
        bottomRight = max(objFeatures.PixelList);
        objFeatures.BoundingBox = [topLeft, bottomRight-topLeft+1];
    end
end

function [major,minor,eccentricity,orientation] = computeEllipseParameters(list,centroid)
% COMPUTEELLIPSEPARAMETERS  Compute parameters of ellipse for object
%
% Taken directly from matlab's regionprops

% Assign X and Y variables so that we're measuring orientation
% counterclockwise from the horizontal axis.

xbar = centroid(1);
ybar = centroid(2);

x = list(:,1) - xbar;
y = -(list(:,2) - ybar); % This is negative for the
% orientation calculation (measured in the
% counter-clockwise direction).

N = length(x);

% Calculate normalized second central moments for the region. 1/12 is
% the normalized second central moment of a pixel with unit length.
uxx = sum(x.^2)/N + 1/12;
uyy = sum(y.^2)/N + 1/12;
uxy = sum(x.*y)/N;

% Calculate major axis length, minor axis length, and eccentricity.
common = sqrt((uxx - uyy)^2 + 4*uxy^2);
major = 2*sqrt(2)*sqrt(uxx + uyy + common);
minor = 2*sqrt(2)*sqrt(uxx + uyy - common);
eccentricity = 2*sqrt((major/2)^2 - (minor/2)^2) / major;

% Calculate orientation.
if (uyy > uxx)
    num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
    den = 2*uxy;
else
    num = 2*uxy;
    den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
end
if (num == 0) && (den == 0)
    orientation = 0;
else
    orientation = (180/pi) * atan(num/den);
end

end

function perimeter = computePerimeterFromBoundary(B)
% COMPUTEPERIMETERFROMBOUNDARY  Compute perimeter
%
% Directly taken from matlab's regionprops

delta = diff(B).^2;
if(size(delta,1) > 1)
    isCorner  = any(diff([delta;delta(1,:)]),2); % Count corners.
    isEven    = any(~delta,2);
    perimeter = sum(isEven)*0.980 + sum(~isEven)*1.406 - sum(isCorner)*0.091;
else
    perimeter = 0; % if the number of pixels is 1 or less.
end
end

function [meanKpos, meanKneg] = computeBoundaryCurvatureMeans(B, options)
% COMPUTEBOUNDARYCURVATUREMEANS  Compute the mean positive and mean
% negative boundary curvature

% James Kapaldo
% 2016-11-18

G = options.G;
dG = options.dG;
ddG = options.ddG;

cntr = mean(B,1);

B = bsxfun(@minus,B,cntr);

% Get boundary curvature
Bs = imfilter(B,G,'circular','conv','same');
d_B = imfilter(Bs,dG,'circular','conv','same');
dd_B = imfilter(Bs,ddG,'circular','conv','same');

d_B(end,:) = [];
dd_B(end,:) = [];

kappa = (d_B(:,1).*dd_B(:,2) - d_B(:,2).*dd_B(:,1))./(sum(d_B.^2,2).^(3/2));
meanKpos = mean(kappa(kappa>0));
meanKneg = mean(kappa(kappa<0));

if isnan(meanKpos)
    meanKpos = 0;
end
if isnan(meanKneg)
    meanKneg = 0;
end

end




















% 
% 
% 
% function out = extractObjectFeatures(CC,I,offset,imSize,featureList,options)
% % Take in connected component structure returned from bwconncomp() and
% % compute several features
% % - area
% % - meanIntensity
% % - centroid
% %
% % should be faster than regionprops if only measuring these
% 
% % 2016-11-16 : added in returnPixelList flag
% % 2016-11-17 : removed returnPixelList flag -- just check for two output
% % arguments
% 
% Nb = options.BoundarySize; % 100
% Nfd = options.NumberFourierDescriptors; % 16
% 
% pixIdx = CC.PixelIdxList;
% N = CC.NumObjects;
% L = labelmatrix(CC);
% 
% % N = max( 0, floor(double(max(L(:)))) );
% % pixIdx = label2idxmex(L, double(N));
% 
% [m,n,ch] = size(I); %#ok<ASGLU>
% 
% 
% % Define some flags
% 
% computeShape = false;
% computeIntensity = false;
% computeRadialIntensity = false;
% computeTexture = false;
% computeGranularity = false;
% computeLBP = false;
% computeCorrelation = false;
% 
% 
% 
% 
% returnSimplifiedBoundary = false;
% returnLinPixelList = false;
% returnBoundingBox = false;
% 
% computePixelList = true;
% 
% numFeatures = 0;
% 
% for i = 1:numel(featureList)
%     switch featureList{i}
%         case 'Shape'
%             % Centroid
%             % Area
%             % Ellipse parameters : 
%             %   major axis length
%             %   minor axis length
%             %   eccentricity
%             %   orientation
%             % Curvature : 
%             %   mean positive curvature
%             %   mean negative curvature
%             % Solidity : area/convexHullArea
%             % IsoperimetricRatio : 4*pi*Area/perimeter^2  -- ratio of the area to the area of a circle with the same perimeter, this is called form factor in cell profiler
%             % FourierDescriptors
%             
%             % Get region boundaries
%             B = bwboundariesmex(double(L), 8); % regionboundariesmex(double(L),8); 
% 
%             % I dont know what the difference between these two functions
%             % are, but I want my boundaries to be clockwise oriented and I
%             % used bwboundaries before with the fourier descriptors
%             % function without having to check the orientation; so, i will
%             % use that mex function.
% 
%             % Create filters for getting derivatives
%             kappaSmoothingSigma = 2;
% 
%             filtSize = round(7*kappaSmoothingSigma);
%             l = -floor(filtSize/2):ceil(filtSize/2);
%             G = exp(-(l).^2/(2*kappaSmoothingSigma^2))/(kappaSmoothingSigma*sqrt(2*pi));
%             dG = -(l) .* G / kappaSmoothingSigma^2;
%             ddG = - G / kappaSmoothingSigma^2 + (l).^2 .* G / kappaSmoothingSigma^4;
%             
%             G = G(:);
%             dG = dG(:);
%             ddG = ddG(:);
%             
%             computeShape = true;
%             
%             shapeFeatureNames = {'Centroid_x',...
%                                  'Centroid_y',...
%                                  'Ellipse_MajorAxisLength',...
%                                  'Ellispe_MinorAxisLength',...
%                                  'Ellipse_Eccentricity',...
%                                  'Ellispe_Orientation',...
%                                  'Curvature_MeanPositive',...
%                                  'Curvature_MeanNegative',...
%                                  'Solidity',...
%                                  'FormFactor',...
%                                  cellstr('FourierDescriptor_'+string(1:Nfd))};
%                              
%             shapeFeatures = zeros(N,11);
%             computePixelList = true;
%             
%         case 'Intensity' % min, mean, max, 1st quartile, median, 3rd quartile, std, mad
%             computeIntensity = true;
%             
%             intensityFeatureNames = ('Intensity_' + string({'Min_', 'Mean_', 'Max_', '1stQuartile_','Median_','3rdQuartile_','Std_','MAD_'}) + string(1:ch)')';
%             intensityFeatureNames = intensityFeatureNames(:)';
%             
%             intensityFeatures = zeros(N,8*ch);
%             
%         case 'RadialIntensity'
%             computeRadialIntensity = true;
%             
%         case 'Texture'
%             computeTexture = true;
%             
%         case 'Granularity'
%             computeGranularity = true;
%             
%         case 'LBP'
%             computeLBP = true;
%             
%         case 'Correlation'
%             computeCorrelation = true;
%             
%         case 'SimplifiedBoundary'
%             returnSimplifiedBoundary = true;
%             
%         case 'pixelList'
%             returnLinPixelList = true;
%             
%         case 'BoundingBox'
%             returnBoundingBox = true;
%             boundingBox = zeros(N,4);
%     end
% end
% 
% 
% if returnLinPixelList
%     pixelList = cell(numel(N),1);
% end
% 
% 
% % Get the intensities of each channel for each object
% s = [];
% for i = 1:ch
%     s = [s,sprintf('I(x+%d*n*m),',i-1)]; %#ok<AGROW>
% end
% s(end) = [];
% s = ['I = cellfun(@(x) cat(3,' s, '),pixIdx,''UniformOutput'',0);'];
% eval(s);
% 
% 
% % I2 = cellfun(@(x) I(x),I(x+,pixIdx,'UniformOutput',0);
% imSize = imSize(1);
% 
% % Compute fourier descriptors --------------------------------------------
% alreadyComputedSimplifiedBoundaries = false;
% if computeShape
%     if returnSimplifiedBoundary
%         [FD, SimplifiedBoundary] = fourierDescriptors(B,NB,Nfd);
%         alreadyComputedSimplifiedBoundaries = true;
%     else
%         FD = fourierDescriptors(B,NB,Nfd);
%     end
% end
% 
% if returnSimplifiedBoundary && ~alreadyComputedSimplifiedBoundaries
%     [~, SimplifiedBoundary] = fourierDescriptors(B,NB,Nfd);
% end
% 
% 
% % Compute features -------------------------------------------------------
% parfor i = 1:N
%     
%     if computePixelList
%         % Compute pixel list
%         y = rem(pixIdx{i}-1, m) + 1;
%         x = (pixIdx{i} - y)/m + 1;
%     end
%         
%     if computeShape
%         
%         % Compute centroid
%         xbar = mean(x);
%         ybar = mean(y);
% 
%         % Compute ellipse parameters
%         [mj,mn,ecc,ori] = computeEllipseParameters(x,y,xbar,ybar);
%         
%         % Compute boundary curvature means
%         [meanKpos, meanKneg] = computeBoundaryCurvatureMeans(B{i}, G, dG, ddG);
% 
%         % Compute area
%         A = numel(pixIdx{i});
%         
%         % Compute convex hull area
%         K = convhull(B{i}(:,1),B{i}(:,2));
%         hull = B{i}(K,:);
%         
%         topLeft = [min(x),min(y)] - 1;
%         r = hull(:,2) - topLeft(2);
%         c = hull(:,1) - topLeft(1);
%         
%         bottomRight = [max(x),max(y)];
%         boundingBox(i,:) = [topLeft+1, (bottomRight-topLeft)];
%         
%         M = bottomRight(2)-topLeft(2);
%         N = bottomRight(1)-topLeft(1);
%         
%         convexImage = roipoly(M,N,c,r);
%         CHA = sum(convexImage(:));
%         
%         % Compute perimeter
%         P = computePerimeterFromBoundary(B{i});
% 
%         % Save the features
%         shapeFeatures(i,:) = [xbar,ybar,A,mj,mn,ecc,ori,meanKpos,meanKneg,A/CHA,4*pi*A/P^2]
%     end
%     
%     if computeIntensity
%         % min, mean, max, 1st quartile, median, 3rd quartile, std, mad
%         minI = min(I{i},[],1); % 1x1x3
%         meanI = mean(I{i},1); % 1x1x3
%         maxI = max(I{i},[],1); % 1x1x3
%         prctileI = permute(prctile(I{i},[25,50,75]),[2,1,3]); % 1x3x3
%         stdI = std(I{i},[],1); % 1x1x3
%         madI = median(I{i}-prctileI(1,2,:),1); % mad(I{i},1); % 1x1x3
%         
%         features = [minI, meanI, maxI, prctileI, stdI, madI]; % should be 1 x L x ch
%         features = permute(features,[2,3,1]);
%         
%         intensityFeatures(i,:) = features(:)';
%     end
%     
%     if returnLinPixelList
%         pixelList{i} = (y+offset(2)) + (x+offset(1)-1)*imSize; %#ok<PFBNS>
%     end
% end
% 
% if computeShape
%     shapeFeatures(:,1:2) = shapeFeatures(:,1:2) + offset;
%     shapeFeatures = [shapeFeatures, FD];
% else
%     shapeFeatures = [];
% end
% 
% if ~computeIntensity
%     intensityFeatures = [];
% end
% 
% features = [shapeFeatures, intensityFeatures];
% boundingBox(:,1:2) = boundingBox(:,1:2) + offset;
% 
% out.features = features;
% if returnLinPixelList
%     out.pixelList = pixelList;
% end
% if returnBoundingBox
%     out.boundingBox = boundingBox;
% end
% if returnSimplifiedBoundary
%     out.SimplifiedBoundary = SimplifiedBoundary;
% end
% 
% end
% 
% 
% 
% function [major,minor,eccentricity,orientation] = computeEllipseParameters(x,y,xbar,ybar)
% % COMPUTEELLIPSEPARAMETERS  Compute parameters of ellipse for object
% %
% % Taken directly from matlab's regionprops
% 
% % Assign X and Y variables so that we're measuring orientation
% % counterclockwise from the horizontal axis.
% 
% % xbar = centroid(1);
% % ybar = centroid(2);
% 
% % x = list(:,1) - xbar;
% % y = -(list(:,2) - ybar); % This is negative for the
% % orientation calculation (measured in the
% % counter-clockwise direction).
% 
% x = x - xbar;
% y = -(y-ybar);
% 
% N = length(x);
% 
% % Calculate normalized second central moments for the region. 1/12 is
% % the normalized second central moment of a pixel with unit length.
% uxx = sum(x.^2)/N + 1/12;
% uyy = sum(y.^2)/N + 1/12;
% uxy = sum(x.*y)/N;
% 
% % Calculate major axis length, minor axis length, and eccentricity.
% common = sqrt((uxx - uyy)^2 + 4*uxy^2);
% major = 2*sqrt(2)*sqrt(uxx + uyy + common);
% minor = 2*sqrt(2)*sqrt(uxx + uyy - common);
% eccentricity = 2*sqrt((major/2)^2 - (minor/2)^2) / major;
% 
% % Calculate orientation.
% if (uyy > uxx)
%     num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
%     den = 2*uxy;
% else
%     num = 2*uxy;
%     den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
% end
% if (num == 0) && (den == 0)
%     orientation = 0;
% else
%     orientation = (180/pi) * atan(num/den);
% end
% 
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function perimeter = computePerimeterFromBoundary(B)
% % COMPUTEPERIMETERFROMBOUNDARY  Compute perimeter
% %
% % Directly taken from matlab's regionprops
% 
% delta = diff(B).^2;
% if(size(delta,1) > 1)
%     isCorner  = any(diff([delta;delta(1,:)]),2); % Count corners.
%     isEven    = any(~delta,2);
%     perimeter = sum(isEven)*0.980 + sum(~isEven)*1.406 - sum(isCorner)*0.091;
% else
%     perimeter = 0; % if the number of pixels is 1 or less.
% end
% end
% 
% function [meanKpos, meanKneg] = computeBoundaryCurvatureMeans(B, G, dG, ddG)
% % COMPUTEBOUNDARYCURVATUREMEANS  Compute the mean positive and mean
% % negative boundary curvature
% 
% % James Kapaldo
% % 2016-11-18
% 
% cntr = mean(B,1);
% 
% B = bsxfun(@minus,B,cntr);
% 
% % Get boundary curvature
% Bs = imfilter(B,G,'circular','conv','same');
% d_B = imfilter(Bs,dG,'circular','conv','same');
% dd_B = imfilter(Bs,ddG,'circular','conv','same');
% 
% d_B(end,:) = [];
% dd_B(end,:) = [];
% 
% kappa = (d_B(:,1).*dd_B(:,2) - d_B(:,2).*dd_B(:,1))./(sum(d_B.^2,2).^(3/2));
% meanKpos = mean(kappa(kappa>0));
% meanKneg = mean(kappa(kappa<0));
% 
% if isnan(meanKpos)
%     meanKpos = 0;
% end
% if isnan(meanKneg)
%     meanKneg = 0;
% end
% 
% end