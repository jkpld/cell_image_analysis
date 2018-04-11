function [userClassified, imagesUsed] = viewNucleiImages(nucleiImages, channelNames, addClassification,addPhase)
% 

% NOTE, you can click on the images to select a pair nuclei. A pair
% nuclei may be a "sister" nuclei of a late telo/early G1 nulcei, or
% the parent nuclei of a micro-nuclei.

userClassified = [];


if isempty(nucleiImages)
    fprintf('No images given!\n')
    return;
end

if nargin < 4
    addPhase = false;
end

if nargin < 3
    addClassification = false;
end


tileSize = 0.75; % [in] , size of nuclei images.

grt = groot;
if addClassification
    numMaxImagesCanShow = floor(sum(grt.MonitorPositions(:,3))/grt.ScreenPixelsPerInch/tileSize - 5);
else
    numMaxImagesCanShow = floor(sum(grt.MonitorPositions(:,3))/grt.ScreenPixelsPerInch/tileSize - 2);
end

K = size(nucleiImages,4);
selIdx = 1:K;
if K > numMaxImagesCanShow
    warning('Only %d of the nuclei images can be shown at a time. I have selected a random %d nuclei out of the given %d nuclei.', numMaxImagesCanShow,numMaxImagesCanShow,numel(selIdx))
    imagesUsed = sort(randperm(K,numMaxImagesCanShow),'ascend');
    nucleiImages = nucleiImages(:,:,:,imagesUsed);
    selIdx = selIdx(imagesUsed);
    K = numMaxImagesCanShow;
else
    imagesUsed = 1:K;
end

if addClassification
%     classificationGroup = {'ringLike', 'condenced','splitting','justSplit','interphase','microNuclei','apoptotic','damage2','blurry','manyNuclei'};
    phaseDamageGroup = {'ringLike', 'condensed','splitting','justSplit','interphase','apoptotic'};
    blurryGroup = {'notBlurry','Blurry'};
    nucleiCountGroup = {'goodSegment','badSegment'};
    clsfctnGroupNames = {'phaseDamage','blurry','segmentation'}; % might also add in subSphase later.
    
    radioGroups = {phaseDamageGroup,blurryGroup,nucleiCountGroup};
    nonRadioGroups = {};
    
    numOptions = 0;
    for i = 1:length(radioGroups)
        numOptions = numOptions + length(radioGroups{i});
    end
    for i = 1:length(nonRadioGroups)
        numOptions = numOptions + length(nonRadioGroups{i});
    end
    
    
    bM = ((10*addPhase + 1 + numOptions)*0.5 + (1+length(radioGroups)+length(nonRadioGroups))*0.1 + 1)/2.54;
    lM = 5/2.54;
else
    bM = 0.2;
    lM = 0.5;
end


maxNumChannels = sum(channelNames~="Mask");

wAx = length(selIdx)*tileSize; % 0.75 inches per nuclei image with 96 pixels per inch
hAx = tileSize;

tM = 0.2;%0.9;
rM = 0.7;

W = wAx + lM + rM;
H = hAx*maxNumChannels + tM + bM;

idx = 1;
pos = [grt.MonitorPositions(idx,1)/grt.ScreenPixelsPerInch + 0.1, ...
    grt.MonitorPositions(idx,4)/grt.ScreenPixelsPerInch - 1.5 - H, ...
    W,H];

viewNucFig = figure('Units','inch','Position',pos,'MenuBar','none','Visible','off');
try
for i = maxNumChannels:-1:1
    pos = [lM, bM+(i-1)*hAx, wAx, hAx];
    ax(i) = axes('Parent',viewNucFig,'Units','inch','Position',pos);
end

%% Create brighten/darken buttons
Icon = [...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0;...
    0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0;...
    0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

Icon = cat(3,Icon,Icon,Icon);
Icon = uint8(Icon*(210-40) + 40);

for i = 2:-1:1
    contrastControl(i).br = uicontrol(viewNucFig,'style','pushbutton',...
        'units','inch',...
        'position',[sum(ax(i).Position([1 3]))+0.4, ax(i).Position(2)+ax(i).Position(4)/2 + 0.06, 16/96, 16/96],...
        'CData',Icon,...
        'UserData',i);
    
    contrastControl(i).dr = uicontrol(viewNucFig,'style','pushbutton',...
        'units','inch',...
        'position',[sum(ax(i).Position([1 3]))+0.4, ax(i).Position(2)+ax(i).Position(4)/2 - 0.04-16/96, 16/96, 16/96],...
        'CData',flipud(Icon),...
        'UserData',i);
end

%% Plot nuclei and information
colors = repmat({[1 1 1]},1,maxNumChannels);

dapiInd = strcmp('DAPI',channelNames);
gfpInd = strcmp('GFP',channelNames);
colors{dapiInd} = [0 0 1];
colors{gfpInd} = [0 1 0];
% colors{dapiInd} = [1 1 1];
% colors{gfpInd} = [1 1 1];
if maxNumChannels == 3
    colors{~(dapiInd | gfpInd)} = [1 0 0];
end


plotNucleiMontage(nucleiImages,channelNames,colors,ax,contrastControl)

%% Plot classification buttons

if addClassification
    classifyAllParent = uicontainer('Parent',viewNucFig,'Units','Centimeter','Position',[0.5 0.5 4,bM*2.54-0.9]);
    [classifyAllButtons, classificationNames, classificationIndex] = createClassificationOptions(classifyAllParent,radioGroups,true, @masterSet,addPhase);
    
    
    for i = K:-1:1
        tmpParent = uicontainer('Parent',viewNucFig,'Units','Centimeter','Position',[lM*2.54+(i-1)*0.75*2.54+0.3, 0.5, 0.5*2.54,bM*2.54-0.9]);
        nucleiClassification(:,i) = createClassificationOptions(tmpParent,radioGroups,false,@slaveSet,addPhase);
    end
    
    I = 0.4*ones(23,79);
    Ic = insertText(I,flip(ceil((size(I)/2)))+1,'Save & Exit','FontSize',13,'BoxColor','black','TextColor','white','BoxOpacity',0.4,'AnchorPoint','center');
    uicontrol('Parent',viewNucFig,'Style','Pushbutton','Units','centimeter',...
        'Position',[1.25,bM*2.54,size(I,2)/96*2.54,size(I,1)/96*2.54],'Callback',@saveExit,'CData',Ic);
    
    uicontrol('Parent',viewNucFig,'Style','text','Units','centimeter',...
        'Position',[0.5,H*2.54-3.5,3.5,3],'String','Classify the nuclei with the buttons below and click Save&Exit when done. If a nuclei has already been clasified, its classifications will show up and any changes will override the previous classification.',...
        'ForegroundColor',0.6*[1 1 1],'HorizontalAlignment','left',...
        'HandleVisibility','off');
    
    
%     [prviouslySavedClsfctns, prviouslySelectedSisterIndex] = getSavedClsfctns(obj,selIdx,clsfctnGroupNames,classificationIndex);
%     for l = 1:numel(prviouslySavedClsfctns)
%         nucleiClassification(l).Value = prviouslySavedClsfctns(l);
%     end
    
    slaveSet;
end

%% Show figure

setTheme(viewNucFig,'dark')
set(ax,'TickDir','out','TickLength',[0.01*5/length(selIdx), 0.025])
viewNucFig.Visible = 'on';

catch ME
    close(viewNucFig);
    rethrow(ME)
end
if nargout > 0 && (addClassification || addPhase)
uiwait;
end
% if nargout > 0
%     varargout{1} = viewNucFig;
% end
%% Create callback functions

    function masterSet(src,group)
%         src
        %         fprintf('masterSet\n')
        if src == 1
            for rb = 1:length(classifyAllButtons)
                set(nucleiClassification(rb,:),'Value',classifyAllButtons(rb).Value)
            end
        else
            set(nucleiClassification(1,:),'Value',classifyAllButtons(1).Value)
            for rb = 1:length(group)
                set(nucleiClassification(group(rb),:),'Value',classifyAllButtons(group(rb)).Value)
            end
        end
            
    end

    function slaveSet(~,~)
        %         fprintf('slaveSet\n')
        for rb = 1:length(classifyAllButtons)
            if ~all([nucleiClassification(rb,:).Value])
                classifyAllButtons(rb).Value = 0;
            else
                classifyAllButtons(rb).Value = 1;
            end
        end
        clsfctn = reshape([nucleiClassification.Value],size(nucleiClassification));
        for k = 1:length(selIdx)
            nucleiClassification(1,k).Value = all(~clsfctn(2:end,k),1);
        end
    end

    function saveExit(~,~)
        
        % This function is a bit slopy and can be made better by only
        % updated the classifications of the nuclei whose classifications
        % actually changed. Determine if the classifications changed by
        % comparing the columns of clsfctn and previouslySavedClsfctns. 
        %
        % Leaving it as is should not waist much time though since we only
        % classify maybe 45 at a time.
        
        clsfctn = reshape([nucleiClassification.Value],size(nucleiClassification));
%         classNames = string([classificationNames{2:end}]);
        
%         clsfctn(1,:) = []; % remove dont't classify
        
        userClassified = struct();
        
        for groupNum = 1:numel(clsfctnGroupNames)
            userClassified.(clsfctnGroupNames{groupNum}) = struct();
            userClassified.(clsfctnGroupNames{groupNum}).name = lower(classificationNames{groupNum+1});
            userClassified.(clsfctnGroupNames{groupNum}).classification = clsfctn(classificationIndex{groupNum+1},:).';
        end

        close(viewNucFig);
    end


end


%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
