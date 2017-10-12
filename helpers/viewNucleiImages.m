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
uiwait;
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
%                 'name',classificationNames{groupNum+1}, ...
%                 'classification',clsfctn(classificationIndex{groupNum+1},:).');
        end
        
%         userClassified.segmentation.classification
%         if ~(isequal(clsfctn, prviouslySavedClsfctns) && isequal(prviouslySelectedSisterIndex,pairNucleiIdx)) % don't need to do anything if nothing changed
%                         
%             dontClsfy = logical(clsfctn(1,:));
%             clsfctnIdx = selIdx(~dontClsfy);
%             sisterIdx = pairNucleiIdx(~dontClsfy);
%             clsfctn = clsfctn(:,~dontClsfy)';
%             dontClsfy = selIdx(dontClsfy);
% %             sisterIdx(sisterIdx==0) = uint8(intmax('uint8'));
%             for groupNum = 1:numel(clsfctnGroupNames)
%                 
%                 groupName = clsfctnGroupNames{groupNum};
%                 groupIndices = classificationIndex{groupNum+1};
%                 
%                 % If the group does not exist, create it.
%                 if ~isfield(userClassified,groupName)
%                     userClassified.(groupName) = struct('index',intmax('uint32')*ones(1000,1,'uint32'),'classification',intmax('uint8')*ones(1000,numel(groupIndices),'uint8'),'name',{{}},'counter',1);
%                 end
%                 
%                 % If the group was just created, add in the classification
%                 % names of the group. If the group was previously created,
%                 % then make sure the new classification names match the old
%                 % classification names
%                 if isempty(userClassified.(groupName).name)
%                     userClassified.(groupName).name = lower(classificationNames{groupNum+1});
%                 else
%                     if ~isequal(lower(userClassified.(groupName).name), lower(classificationNames{groupNum+1}))
%                         error('veiwTrain:badClassificationNames', 'The classifications names in group %s do not match the names already in that group.', groupName)
%                     end
%                 end
%                 
%                 % Get the classifications for this group
%                 groupClsfctn = clsfctn(:,groupIndices); % numNuclei x classificationOptions
%                 
%                 % remove the nuclei not classified
%                 notClsfd = all(groupClsfctn==0,2);
%                 clsfdNucleiIndices = uint32(clsfctnIdx(~notClsfd));
%                 clsfdSisterIndices = uint32(sisterIdx(~notClsfd));
%                 numClsfdNuclei = numel(clsfdNucleiIndices);
%                 groupClsfctn(notClsfd,:) = [];
%                 
%                 % get the next index of the next element to be saved.
%                 counter = userClassified.(groupName).counter;
%                 
%                 % add to the size of the index and value arrays if
%                 % necessary  -- i know this may add the extra elements a
%                 % bit before they are actually needed, but it should be
%                 % fine
%                 if (numClsfdNuclei + counter) > numel(userClassified.(groupName).index)
%                     
%                     userClassified.(groupName).index = [userClassified.(groupName).index; intmax('uint32')*ones(1000,1,'uint32')];
%                     userClassified.(groupName).classification = [userClassified.(groupName).classification; intmax('uint8') *ones(1000,numel(groupIndices),'uint8')];
%                     
%                     if strcmp(groupName,'phase')
%                         userClassified.(groupName).sister_index = [userClassified.(groupName).sister_index; intmax('uint32')*ones(1000,1,'uint32')];
%                     end
%                 end
%                 
%                 % before we save the new values, we need check if these
%                 % nuclei have already been assigned. if they have been
%                 % assigned and have the same value, fine; however, if they
%                 % have been saved and have a different value, then we need
%                 % to throw up a message asking the user to pick one.
%                 % -- actually, just overwrite the previous classificaiton
%                 
%                 for j = 1:numClsfdNuclei
%                     
%                     prviousClsfctn = find(clsfdNucleiIndices(j) == userClassified.(groupName).index(1:counter-1));
%                     
%                     if ~isempty(prviousClsfctn)
%                         
%                         userClassified.(groupName).classification(prviousClsfctn,:) = groupClsfctn(j,:);
%                         
%                         if strcmp(groupName,'phase')
%                             userClassified.(groupName).sister_index(prviousClsfctn) = clsfdSisterIndices(j);
%                         end
%                     else
%                         userClassified.(groupName).index(counter) = clsfdNucleiIndices(j);
%                         userClassified.(groupName).classification(counter,:) = groupClsfctn(j,:);
%                         if strcmp(groupName,'phase')
%                             userClassified.(groupName).sister_index(counter) = clsfdSisterIndices(j);
%                         end
%                         counter = counter + 1;
%                     end
%                 end
% 
%                 for j = 1:numel(dontClsfy)
%                     
%                     prviousClsfctn = find(dontClsfy(j) == userClassified.(groupName).index(1:counter-1));
%                     
%                     if ~isempty(prviousClsfctn)
%                         
%                         userClassified.(groupName).index(prviousClsfctn) = [];
%                         userClassified.(groupName).classification(prviousClsfctn,:) = [];
%                         if strcmp(groupName,'phase')
%                             userClassified.(groupName).sister_index(prviousClsfctn) = [];
%                         end
%                         counter = counter - 1;
%                     end
%                 end
%                 
%                 userClassified.(groupName).counter = counter;
%             end
%         end
        close(viewNucFig);
    end


end


% 
% function saveExit(~,~)
%         
%         % This function is a bit slopy and can be made better by only
%         % updated the classifications of the nuclei whose classifications
%         % actually changed. Determine if the classifications changed by
%         % comparing the columns of clsfctn and previouslySavedClsfctns. 
%         %
%         % Leaving it as is should not waist much time though since we only
%         % classify maybe 45 at a time.
%         
%         clsfctn = reshape([nucleiClassification.Value],size(nucleiClassification));
%         
%         if ~isequal(clsfctn, prviouslySavedClsfctns)
%                         
%             dontClsfy = logical(clsfctn(1,:));
%             clsfctnIdx = selIdx(~dontClsfy);
%             sisterIdx = pairNucleiIdx(~dontClsfy);
%             clsfctn = clsfctn(:,~dontClsfy)';
%             dontClsfy = selIdx(dontClsfy);
%             
%             for groupNum = 1:numel(clsfctnGroupNames)
%                 
%                 groupName = clsfctnGroupNames{groupNum};
%                 
%                 % If the group does not exist, create it.
%                 if ~isfield(obj.userClassified,groupName)
%                     obj.userClassified.(groupName) = struct('index',intmax('uint32')*ones(1000,1,'uint32'),'value',intmax('uint8')*ones(1000,1,'uint8'),'name',{{}},'counter',1);
%                     if strcmp(groupName,'phase')
%                         obj.userClassified.(groupName).sister_index = intmax('uint32')*ones(1000,1,'uint32');
%                     end
%                 end
%                 
%                 % If the group was just created, add in the classification
%                 % names of the group. If the group was previously created,
%                 % then make sure the new classification names match the old
%                 % classification names
%                 if isempty(obj.userClassified.(groupName).name)
%                     obj.userClassified.(groupName).name = lower(classificationNames{groupNum+1});
%                 else
%                     if ~isequal(lower(obj.userClassified.(groupName).name), lower(classificationNames{groupNum+1}))
%                         error('veiwTrain:badClassificationNames', 'The classifications names in group %s do not match the names already in that group.', groupName)
%                     end
%                 end
%                 
%                 
%                 % Get the classifications for this group
%                 groupClsfctn = clsfctn(:,classificationIndex{groupNum+1});
%                 
%                 % remove the nuclei not classified
%                 notClsfd = all(groupClsfctn==0,2);
%                 groupClsfctnIdx = uint32(clsfctnIdx(~notClsfd));
%                 groupSisterIdx = uint32(sisterIdx(~notClsfd));
%                 numClsfctns = numel(groupClsfctnIdx);
%                 groupClsfctn(notClsfd,:) = [];
%                 
%                 % get the next index of the next element to be saved.
%                 counter = obj.userClassified.(groupName).counter;
%                 
%                 % add to the size of the index and value arrays if
%                 % necessary  -- i know this may add the extra elements a
%                 % bit before they are actually needed, but it should be
%                 % fine
%                 if (numClsfctns + counter) > numel(obj.userClassified.(groupName).index)
%                     obj.userClassified.(groupName).index = [obj.userClassified.(groupName).index; intmax('uint32')*ones(1000,1,'uint32')];
%                     obj.userClassified.(groupName).classification = [obj.userClassified.(groupName).classification; intmax('uint8') *ones(1000,1,'uint8')];
%                     if strcmp(groupName,'phase')
%                         if ~isfield(obj.userClassified.(groupName),'sister_index')
%                             obj.userClassified.(groupName).sister_index = intmax('uint32')*ones(numel(obj.userClassified.(groupName).index),1,'uint32');
%                         else
%                             obj.userClassified.(groupName).sister_index = [obj.userClassified.(groupName).sister_index; intmax('uint32')*ones(1000,1,'uint32')];
%                         end
%                     end
%                 end
%                 
%                 % create the classification values
%                 groupClsfctnValue = uint8(sum(bsxfun(@times,groupClsfctn,0:size(groupClsfctn,2)-1),2));
%                 
%                 
%                 % before we save the new values, we need check if these
%                 % nuclei have already been assigned. if they have been
%                 % assigned and have the same value, fine; however, if they
%                 % have been saved and have a different value, then we need
%                 % to throw up a message asking the user to pick one.
%                 % -- actually, just overwrite the previous classificaiton
%                 
%                 for j = 1:numClsfctns
%                     prviousClsfctn = find(groupClsfctnIdx(j) == obj.userClassified.(groupName).index(1:counter-1));
%                     if ~isempty(prviousClsfctn)
%                         obj.userClassified.(groupName).classification(prviousClsfctn) = groupClsfctnValue(j);
%                         if strcmp(groupName,'phase') && ~isnan(groupSisterIdx(j))
%                             if ~isfield(obj.userClassified.(groupName),'sister_index')
%                                 obj.userClassified.(groupName).sister_index = intmax('uint32')*ones(numel(obj.userClassified.(groupName).index),1,'uint32');
%                             end
%                             obj.userClassified.(groupName).sister_index(prviousClsfctn) = groupSisterIdx(j);
%                         end
%                     else
%                         obj.userClassified.(groupName).index(counter) = groupClsfctnIdx(j);
%                         obj.userClassified.(groupName).classification(counter) = groupClsfctnValue(j);
%                         if strcmp(groupName,'phase') && ~isnan(groupSisterIdx(j))
%                             if ~isfield(obj.userClassified.(groupName),'sister_index')
%                                 obj.userClassified.(groupName).sister_index = intmax('uint32')*ones(numel(obj.userClassified.(groupName).index),1,'uint32');
%                             end
%                             obj.userClassified.(groupName).sister_index(counter) = groupSisterIdx(j);
%                         end
%                         counter = counter + 1;
%                     end
%                 end
% 
%                 for j = 1:numel(dontClsfy)
%                     prviousClsfctn = find(dontClsfy(j) == obj.userClassified.(groupName).index(1:counter-1));
%                     if ~isempty(prviousClsfctn)
%                         obj.userClassified.(groupName).index(prviousClsfctn) = [];
%                         obj.userClassified.(groupName).classification(prviousClsfctn) = [];
%                         if strcmp(groupName,'phase')
%                             if isfield(obj.userClassified.(groupName),'sister_index')
%                                 obj.userClassified.(groupName).sister_index(prviousClsfctn) = [];
%                             end
%                         end
%                         counter = counter - 1;
%                     end
%                 end
%                 
%                 obj.userClassified.(groupName).counter = counter;
%             end
%         end
%         close(viewNucFig);
%     end