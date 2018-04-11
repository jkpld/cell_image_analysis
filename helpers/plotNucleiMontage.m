function plotNucleiMontage(images,channelNames,colors,ax,contrastControl)
% images  -  an image array with the same format as a grayscale image
%            accepted by the function montage accepts. it is assumed that
%            all images of one channel are next to each onther in the
%            array, i.e. if there are two channels images(:,:,1,1:n/2) are
%            channel 1, and images(:,:,1,n/2+1:n) are channel 2
%
% channelsNames  -  cell array whose length is equal to the number of
% channels and whose elements are the names of the channels. 
%
% colors  -  cell array with number of elements equal to numChannels that
% gives the colors for each channel. colors can be empty to give black and
% white images
%
% parents  -  array of axis parents with the number of elements equal to
% numChannels

if nargin < 5
    contrastControl = [];
end

imgMin = min(min(images,[],1),[],2);
imgMax = max(max(images,[],1),[],2);

images = images - imgMin;
images = images./(imgMax-imgMin);

[Y,X,C,K] = size(images);
images = reshape(permute(images,[1,2,4,3]),[Y,X*K,C]);

mIdx = channelNames=="Mask";
if any(mIdx)
    mask = images(:,:,mIdx);
    CC = bwconncomp(mask);
    objCenter = regionprops(CC,'Centroid');
    objCenter = {objCenter.Centroid};
    objCenter = cat(1,objCenter{:});
    [~,objind] = pdist2(objCenter, [(ceil(X/2):X:X*K)',ceil(Y/2)*ones(K,1)],'euclidean','Smallest',1);

    objMaskCC = CC.PixelIdxList(objind);
    objMask = zeros(Y,X*K,'like',mask);
    objMask(cat(1,objMaskCC{:})) = true;
    allB = bwboundaries(mask);
    oB = bwboundaries(objMask);

    allB = cellfun(@(x) [x;nan,nan],allB,'UniformOutput',0);
    oB = cellfun(@(x) [x;nan,nan],oB,'UniformOutput',0);
    allB = cat(1,allB{:});
    oB = cat(1,oB{:});
    C = C-1;
end

images = images(:,:,~mIdx);
images = reshape(permute(images,[2,1,3]),[X*K,Y*C]).';

lims = cell(1,C);
brightnessLevels(C) = 0;

for i = 1:C
    I = round(single(intmax('uint16'))*images((1:Y) + (i-1)*Y, :));
    
    ax(i).XLim = 0.5 + [0, K*size(images,1)];
    ax(i).XLim = 0.5 + [0, size(images,1)];

    lims{i} = [round(prctile(I(:),80:2:100)),single(intmax('uint16'))];
    brightnessLevels(i) = numel(lims{i})-2;
    
    image(I,'Parent',ax(i),'HitTest','off')

    if ~isempty(colors)
        col = colors{i};
        cmap = generateColorGrad([0 0 0],col,lims{i}(brightnessLevels(i)));
        colormap(ax(i),cmap);
    else
        col = [0,1,1];
    end
    if isequal(colors{i},[1,1,1])
        col = [0,1,1];
    end
    ax(i).YLabel.String = channelNames{i};
    
    c = colorbar(ax(i));
    c.Units = 'inch';
    c.Position = [sum(ax(i).Position([1,3]))+0.02, ax(i).Position(2)+0.05, 0.1, ax(i).Position(4)-0.1];
    drawnow;
    c.Ruler.Axle.Visible = 'off';
%     c.Ruler.SecondaryLabel.Position(2) = c.Ruler.SecondaryLabel.Position(2)/2;
%     c.Ruler.SecondaryLabel.Position(1) = 3;
%     c.Ruler.SecondaryLabel.VerticalAlignment = 'middle';
    c.Ruler.SecondaryLabel.Visible = 'off';
    c.Ruler.Exponent = 4;
    
    if ~isempty(contrastControl)
        contrastControl(i).br.Callback = @brightenIm;
        contrastControl(i).dr.Callback = @darkenIm;
    end
    
    if any(mIdx)
        line(allB(:,2),allB(:,1),'linestyle',':','color',~col,'Parent',ax(i))
        line(oB(:,2),oB(:,1),'linestyle','--','color',~col,'Parent',ax(i))
    end
end

set(ax,'XTick',1:X:(X*K+1),'XTickLabels','','XGrid','on','Layer','top','Box','on','YTick',[])
        
drawnow;
for i = 1:2
    ax(i).XRuler.Axle.Visible = 'off';
    ax(i).YRuler.Axle.Visible = 'off';
end
set(ax,'Layer','top','TickDir','in')
set(ax,'GridColor',[1 1 1]*0.5,'GridAlpha',0.5)

    function brightenIm(obj,~)
%         init = levels(obj.UserData);
        
        brightnessLevels(obj.UserData) = max(1,min(numel(lims{obj.UserData}),brightnessLevels(obj.UserData)-1));
%         if init ~= levels(obj.UserData)
%             fprintf('brighten\n')
%         end
        if ~isempty(colors)
            cmap = generateColorGrad([0 0 0],colors{obj.UserData},lims{obj.UserData}(brightnessLevels(obj.UserData)));
            colormap(ax(obj.UserData),cmap);
        end
        
    end

    function darkenIm(obj,~)
%         init = levels(obj.UserData);

        brightnessLevels(obj.UserData) = min(numel(lims{obj.UserData}),brightnessLevels(obj.UserData)+1);
%         if init ~= levels(obj.UserData)
%             fprintf('darken\n')
%         end
        if ~isempty(colors)
            cmap = generateColorGrad([0 0 0],colors{obj.UserData},lims{obj.UserData}(brightnessLevels(obj.UserData)));
            colormap(ax(obj.UserData),cmap);
        end
    end
end

function c = generateColorGrad(c1,c2,n)
c(n,3) = 0;

for i = 1:3
    c(:,i) = linspace(c1(i),c2(i),n);
end
end
%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
