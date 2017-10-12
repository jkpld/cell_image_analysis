function [buttonGroup, classificationNames,classificationIndex] = createClassificationOptions(parent,radioGroups,addText,extraCallback,addPhase)
%radioGroups is a cell array with each element being a list of names
%from which only one can be selected. The list of names should be a cell
%array.
%nonRadioGroups is a cell array with each element being a list of names
%from which any number may be selected. --- removed
%These groups will be in addition to the standard options 'Dont classify'
%and the cell cycle phases.

% try close(fig5), catch, end;
% fig5 = figure('Position',[1173 551 196 383]);

% addText = true;
if nargin < 5
    addPhase = true;
end
if nargin < 4
    extraCallback = [];
end
nonRadioGroups = {};

names = {'Dont''t classify'};


if addPhase
    phase = {'mitosis', 'ringlike', 'condenced','splitting','justsplit','interphase','g1','s','g2','micronuclei'};
    classificationNames = [names, phase, radioGroups, nonRadioGroups];
    names = [names, phase];
    
    classificationIndex = {1,2:11};
    nxtClsfctnIdx = 11;

    cellCycleLink = [...
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1;...
    2, 0, 1, 1, 1, 1, 1, 1, 1, 1;...
    2, 1, 0, 1, 1, 1, 1, 1, 1, 1;...
    2, 1, 1, 0, 1, 1, 1, 1, 1, 1;...
    2, 1, 1, 1, 0, 1, 1, 1, 1, 1;...
    1, 1, 1, 1, 1, 0, 1, 1, 1, 1;...
    1, 1, 1, 1, 1, 2, 0, 1, 1, 1;...
    1, 1, 1, 1, 1, 2, 1, 0, 1, 1;...
    1, 1, 1, 1, 1, 2, 1, 1, 0, 1;...
    1, 1, 1, 1, 1, 1, 1, 1, 1, 0];

    buttonLink = blkdiag(0,cellCycleLink);
    
    gapInd = [1,10];

else
    classificationNames = [names, radioGroups, nonRadioGroups];
    
    classificationIndex = {1};
    nxtClsfctnIdx = 1;
    gapInd = 1;
    
    buttonLink = blkdiag(0);
end






if ~isempty(radioGroups)
    numRGElements = zeros(1,length(radioGroups));
    for i = 1:length(radioGroups)
        names = [names,radioGroups{i}]; %#ok<AGROW>
        bttnLink = ~diag(ones(length(radioGroups{i}),1));
        buttonLink = blkdiag(buttonLink,bttnLink);
        numRGElements(i) = length(radioGroups{i});
        classificationIndex{end+1} = (1:numRGElements(i)) + nxtClsfctnIdx;
        nxtClsfctnIdx = nxtClsfctnIdx + numRGElements(i);
    end
else
    numRGElements = [];
end

if ~isempty(nonRadioGroups)
    numNRGElements = zeros(1,length(nonRadioGroups));
    for i = 1:length(nonRadioGroups)
        names = [names,nonRadioGroups{i}]; %#ok<AGROW>
        bttnLink = zeros(length(nonRadioGroups{i}));
        buttonLink = blkdiag(buttonLink,bttnLink);
        numNRGElements(i) = length(nonRadioGroups{i});
        classificationIndex{end+1} = (1:numNRGElements(i)) + nxtClsfctnIdx;
        nxtClsfctnIdx = nxtClsfctnIdx + numNRGElements(i);
    end
else
    numNRGElements = [];
end
buttonLink(2:end,1) = 1;
buttonLink(1,2:end) = 1;

gapInd = [gapInd,numRGElements,numNRGElements];
gapInd(end) = [];
gapInd = cumsum(gapInd) + 1;

dy = 0.5;
dx = 0.5;
if addText
    lM = 2.5;
else
    lM = 0.25;
end
ddy = 0.1;
N = length(names);

H = dy*N + (1+length(radioGroups)+length(nonRadioGroups))*ddy;

positions = [H-dy*(1:N)', lM*ones(N,1)];
subNodes = zeros(N,1);
gaps = zeros(N,1);
if addPhase
    subNodes([3:6,8:10]) = 1;
end
gaps(gapInd) = 1;

positions(:,1) = positions(:,1)-cumsum(gaps)*ddy;
positions(:,2) = positions(:,2)+subNodes*dx;

% ANNOTATIONS ARE VERY SLOW TO CREATE! MUCH FASTER TO CREATE AN
% AXIS FOR THE SEPERATION LINES AND TEXT! -- LIKE 1 SECOND FASTER
% IN TOTAL
parentSize = parent.Position(3:4);
posLines = (positions(logical(gaps),1) + dy+0.05)/parentSize(2);
tmpAx = axes('Parent',parent','Position',[0 0 1 1]);
X = nan(1,length(posLines)*3);
Y = X;
for i = 1:length(posLines)
    X((1:2)+(i-1)*3) = [0.05,lM+0.5+1*dx]/parentSize(1);
    Y((1:2)+(i-1)*3) = posLines(i);
%     annotation(parent,'line','X',[0.05,lM+0.5+1*dx]/parentSize(1),'Y',posLines(i)*[1 1],'Color',0.5*[1 1 1])
end

line(X,Y,'LineStyle','-','Color',0.5*[1 1 1],'HandleVisibility','off')
tmpAx.Visible = 'off';
tmpAx.HandleVisibility = 'off';
tmpAx.Layer = 'bottom';
tmpAx.XLim = [0 1];
tmpAx.YLim = [0 1];


for i = N:-1:1
    btStyle = 'radiobutton';
    if ~isempty(numNRGElements)
        if i > N-sum(numNRGElements)
            btStyle = 'checkbox';
        end
    end
    buttonGroup(i) = uicontrol(parent,'Style',btStyle,...
        'String','',...
        'Units','centimeter',...
        'Position',[positions(i,2), positions(i,1), 0.5, 0.5],...
        'HandleVisibility','off',...
        'Callback', @classificationChange,...
        'UserData',i);
    if addText
        
%         annotation(parent,'textbox',...
%             'HorizontalAlignment','left',...
%             'VerticalAlignment','middle',...
%             'Units','centimeter',...
%             'HandleVisibility','off',...
%             'EdgeColor','none',...
%             'Position',[positions(i,2)-2.4,positions(i,1),3,0.5],...
%             'String',names{i});
        text((positions(i,2)-2.4)/parentSize(1),(positions(i,1))/parentSize(2)+0.01,names{i},'Parent',tmpAx,...
            'HorizontalAlignment','left',...
            'VerticalAlignment','bottom',...
            'HandleVisibility','off',...
            'EdgeColor','none');
    end
end
buttonGroup(1).Value = 1;


    function classificationChange(obj,~)
        toTurnOff = buttonLink(obj.UserData,:)==1;
        toTurnOn = buttonLink(obj.UserData,:)==2;
        
        set(buttonGroup(toTurnOff),'Value',0);
        set(buttonGroup(toTurnOn),'Value',1);
        
        if ~isempty(extraCallback)
            feval(extraCallback,obj.UserData,classificationIndex{cellfun(@(x) any(x==obj.UserData),classificationIndex)})
        end
    end

end