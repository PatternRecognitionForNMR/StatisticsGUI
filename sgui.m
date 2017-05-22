
function varargout = sgui(varargin)
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @sgui_OpeningFcn, ...
    'gui_OutputFcn',  @sgui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end




function sgui_OpeningFcn(hObject, ~, handles, varargin)
set(hObject,'windowbuttonmotionfcn',@mousemove);

set(handles.view_data,'visible','off')
set(handles.axes1,'visible','off')
set(handles.axes2,'visible','off')
set(handles.axes3,'visible','off')
set(handles.axes4,'visible','off')
set(handles.axes6,'visible','off')

set(handles.DTWButton,'Visible','off');
set(handles.PCAButton,'Visible','off');

handles.output = hObject;
set(handles.pushbutton2,'Enable','off') ;
set(handles.pushbutton3,'Enable','off') ;
set(handles.pushbutton4,'Enable','off') ;
set(handles.cutDataButton,'Enable','off') ;
set(handles.alignData,'Enable','off');

%Create tab group
handles.tgroup = uitabgroup('Parent', handles.figure1,'TabLocation', 'top');

handles.tab1 = uitab('Parent', handles.tgroup, ...
    'Title', 'Data Handling');

handles.tab2 = uitab('Parent', handles.tgroup, ...
    'Title', 'Setup');

handles.tab3 = uitab('Parent', handles.tgroup, ...
    'Title', 'PCA/DTW');

handles.tab4 = uitab('Parent', handles.tgroup, ...
    'Title', 'LDA');

%Place panels into each tab
set(handles.p1,'Parent',handles.tab1)
set(handles.p2,'Parent',handles.tab2)
set(handles.p3,'Parent',handles.tab3)
set(handles.p4,'Parent',handles.tab4)
%Reposition each panel to same location as panel 1
set(handles.p2,'position',get(handles.p1,'position'));
set(handles.p3,'position',get(handles.p1,'position'));
set(handles.p4,'position',get(handles.p1,'position'));

guidata(hObject, handles);


function varargout = sgui_OutputFcn(~, ~, handles)
varargout{1} = handles.output;


function mousemove(hObject,~)
handles = guidata(hObject);
C = get(handles.axes3,'currentpoint');
xlim = get(handles.axes3,'xlim');
ylim = get(handles.axes3,'ylim');
outX = ~any(diff([xlim(1) C(1,1) xlim(2)])<0);
outY = ~any(diff([ylim(1) C(1,2) ylim(2)])<0);
activeTab = find(handles.tgroup.SelectedTab == [handles.tab1 handles.tab2 handles.tab3 handles.tab4 ]);
if outX&&outY && activeTab == 3
    if isfield(handles, 'gspbrush')
        set(handles.gspbrush,'Enable','on')
    end
else
    if isfield(handles, 'gspbrush')
        set(handles.gspbrush,'Enable','off')
    end
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, handles) %#ok<*DEFNU>
% hObject    handle to pushbutton1 (see GCBO)
[FileName,FilePath] = uigetfile({'*.dat;*.mat;*.txt;*.xlsx'});
set(handles.edit2,'string',FileName);


[data_raw, fnames, alldata, steps] = readStatFile([FilePath FileName],get(handles.checkbox5,'value'));
set(handles.cutXVal, 'string',num2str(min(steps(:))));
set(handles.upLim, 'string',num2str(max(steps(:))));
handles.fnames= fnames;
handles.alldata = alldata;
handles.reddata = alldata;
handles.steps = steps;
handles.redsteps = steps;
guidata(hObject,handles);

% Display the size of data
set(handles.text5,'string',size(alldata,2));
set(handles.text6,'string',size(alldata,1));
set(handles.pushbutton2,'Enable','on');
set(handles.pushbutton3,'Enable','on');
set(handles.cutDataButton,'Enable','on');
set(handles.alignData,'Enable','on');


%% draw data in axis
axes(handles.view_data);
cla reset
displaydata = handles.reddata;
colormap(dlmread('wb.map'));
imagesc(linspace(min(steps(:)),max(steps(:)),size(displaydata,1)),1:1:size(displaydata,2),displaydata')
set(handles.view_data,'visible','on');



function edit2_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Reduce by factor edit3
function pushbutton3_Callback(hObject, ~, handles)
factor = get(handles.edit3,'string');
factor = floor(str2double(factor));

if isempty(factor)
    msgbox('Please select integer factor for dataset reduction');
    return
elseif isscalar(factor)~=1 || factor<0
    msgbox('Only positive integers allowed!');
    return
else
    set(handles.edit3, 'string', factor);
end

%get data from handles
alldata = handles.alldata;
steps = handles.steps;

if factor == 1 || factor == 0
    reddata = alldata;
    redsteps = steps;
else
    rest = mod(size(alldata,1),factor); %check if not cleanly devisable by 2
    for i = 1:rest
        alldata(1,:) = [];    %remove one line
        steps(1,:) = [];      %also from steps
    end
    [reddata,redsteps] = deal(zeros(floor(size(alldata,1)/factor),size(alldata,2)));
    for i = 1:size(alldata,2)  %reduce data by averaging
        reddata(:,i) =mean(reshape(alldata(:,i),factor,length(alldata)/factor))';
        redsteps(:,i) = mean(reshape(steps(:,i),factor,length(steps)/factor))';
    end
end

set(handles.text5,'string',size(reddata,2));
set(handles.text6,'string',size(reddata,1));

axes(handles.view_data);
cla reset
colormap(dlmread('wb.map'));
imagesc(redsteps(:,1),1:1:size(reddata,2),reddata')
set(handles.view_data,'visible','on');

handles.reddata = reddata;
handles.redsteps = redsteps;
guidata(hObject,handles);


% --- cutdatabutton
function cutDataButton_Callback(hObject, ~, handles)
if(isfield(handles, 'redsteps'))
    steps = handles.redsteps;
else
    steps = handles.steps;
end

data = handles.reddata;
lowLim = get(handles.cutXVal,'string');
lowLim = floor(str2double(lowLim));

upLim = get(handles.upLim,'string');
upLim = floor(str2double(upLim));

for i = 1 : size(data,2)
    buffer=data(:,i);
    buffer2 = steps(:,i);
    semicutData(:,i) = buffer(buffer2 >= lowLim);
    seminewsteps(:,i) = buffer2(buffer2>=lowLim);
    
    buffer3 = semicutData(:,i);
    buffer4 = seminewsteps(:,i);
    cutData(:,i) = buffer3(buffer4 <= upLim);
    newsteps(:,i) = buffer4(buffer4 <= upLim);
end

axes(handles.view_data);
cla reset
colormap(dlmread('wb.map'));
imagesc(linspace(min(newsteps(:)),max(newsteps(:)),size(cutData,1)),1:1:size(cutData,2),cutData')
set(handles.view_data,'visible','on');
handles.reddata = cutData;
handles.redsteps = newsteps;

set(handles.text5,'string',size(handles.reddata,2));
set(handles.text6,'string',size(handles.reddata,1));

guidata(hObject,handles);

function edit2_Callback(~, ~, ~)

%  press button reduce when pressing enter.
function edit3_Callback(hObject, ~, handles)
pushbutton3_Callback(hObject,0, handles)


function edit3_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Execute evaluation for brushed data
function pushbutton4_Callback(hObject, ~, handles)

redsteps = handles.redsteps;
reddata = handles.reddata;

xd = get(handles.gsp, 'XData');
yd = get(handles.gsp, 'YData');
br = get(handles.gsp, 'BrushData');

brd = [];
xc = [];
yc = [];

if iscell(br)
    for i  = 1 : size(br,1)
        brd = horzcat(brd, br{i});
        xc = horzcat(xc, xd{i});
        yc = horzcat(yc, yd{i});
    end
else
    for i  = 1 : size(br,1)
        brd = horzcat(brd, br(i,:));
        xc = horzcat(xc, xd(i,:));
        yc = horzcat(yc, yd(i,:));
    end
end

xy = [xc ;yc]';
Xproj = [handles.Xproj(:,str2double(handles.pm1)) handles.Xproj(:,str2double(handles.pm2))];

for i= 1: size(reddata,2)
    asgn_fnames(:,i)=find(xy(:,1)==Xproj(i,1)& xy(:,2)==Xproj(i,2),1);
end

axes(handles.axes4);
cla
plot(redsteps(:,asgn_fnames(brd == 1)),reddata(:,asgn_fnames(brd == 1)));
ylabel('amplitude')
xlabel('profile depth')
set(handles.axes4,'visible','on')
legend(handles.fnames(asgn_fnames(brd==1)),'Location','SouthEast');
handles.asgn_fnames = asgn_fnames;
handles.brd = brd;
guidata(hObject,handles)
set(handles.DTWButton,'Visible','on');
set(handles.PCAButton,'Visible','on');




%% procedure selections

% DTW yes/no
function checkbox1_Callback(~, ~, ~)
% PCA yes/no
function checkbox2_Callback(~, ~, ~)
% LDA yes/no
function checkbox3_Callback(~, ~, ~)
% combined sets DTW PCA LDA on.
function checkbox5_Callback(~, ~, ~)
% Normalize Data

function checkbox4_Callback(~, ~, handles)
%% combine yes/no
set(handles.checkbox1, 'value', 1);
set(handles.checkbox2, 'value', 1);
set(handles.checkbox3, 'value', 1);

%% chose components to plot

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
handles.pm1 = str{val};
pushbutton2_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

function popupmenu1_CreateFcn(hObject, ~, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'value', 1);
set(hObject,'Enable','off')

str = get(hObject, 'String');
val = get(hObject,'Value');
handles.pm1 = str{val};
guidata(hObject,handles)


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
handles.pm2 = str{val};
pushbutton2_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, ~, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'value', 2);
set(hObject,'Enable','off')

str = get(hObject, 'String');
val = get(hObject,'Value');
handles.pm2 = str{val};
guidata(hObject,handles)

%% Compute.
function pushbutton2_Callback(hObject, ~, handles)
swDTW = get(handles.checkbox1,'value');
swPCA = get(handles.checkbox2,'value');
swLDA = get(handles.checkbox3,'value');

if(isfield(handles, 'reddata'))
    reddata = handles.reddata;
    fnames = handles.fnames;
else
    msgbox('Please load data first.');
    return
end



%run DTW if checked (1)
if swDTW == 1
    h = waitbar(0,'Initializing waitbar...');
    dist = zeros(length(fnames),length(fnames));
    pdist= zeros(length(fnames),length(fnames));
    
    
    for i= 1:length(fnames)
        h=waitbar(i/length(fnames),h,'running...');
        for j= 1:length(fnames)
            pdist(i,j) = sum(abs(reddata(:,i)-reddata(:,j)),1);
            [dist(i,j)] = dtw(  reddata(:,i)',...
                reddata(:,j)');
        end
    end
    close(h)
    clear i j
    axes(handles.axes1);
    cla
    histogram(nonzeros(triu(pdist)));
    legend('before DTW')
    
    axes(handles.axes2);
    cla
    histogram(nonzeros(triu(dist)));
    legend('after DTW')
    
    set(handles.axes1,'visible','on')
    set(handles.axes2,'visible','on')
    handles.pdist = pdist;
    guidata(hObject,handles)
end


if swPCA ==1
    set(handles.pcarun,'Visible','on');
    set(handles.popupmenu1,'Enable','on')
    set(handles.popupmenu2,'Enable','on')
    
    classes = zeros(1,length(fnames));
    classes(:) = 1;
    
    
    X = reddata';
    [~ , Wpca] = mpca(X);
    Xm = bsxfun(@minus, X, mean(X));
    Xproj = Xm * Wpca;
    
    axes(handles.axes3);
    cla
    gsp = gscatter(Xproj(:,str2double(handles.pm1)),Xproj(:,str2double(handles.pm2)),classes,'k');
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'color', [0 0 0]);  %x-axis
    line(xL, [0 0],'color', [0 0 0]);  %y-axis
    box on
    set(handles.axes3,'visible','on')
    handles.Xproj = Xproj;
    legend('off')
    
    for i = 1:max(classes)
        set(gsp(i), 'Marker','s','MarkerSize', 6,'LineWidth', 1)
    end
    
    set(handles.pushbutton4,'Enable','on') ;
    gspbrush = brush;
    set(gspbrush,'Color','r','Enable','on');
    
    
    
    handles.gsp = gsp;
    handles.gspbrush = gspbrush;
    handles.gspax = gca;
    set(handles.pcarun,'Visible','off');
    guidata(hObject,handles)
end

if swLDA ==1
    X = handles.reddata';
    f = uifigure;
    names = handles.fnames';
    for i = 1 : length(handles.fnames)
        names{i,2} = 1;
    end
    columnname =   {'Name', 'Class'};
    classes_tab = uitable(f,'Data',names,...
        'ColumnName', columnname,...
        'ColumnEditable', [false  true ],...
        'DeleteFcn','classes = TabDeleteFcn(gcbo)');
    
    waitfor(classes_tab)
    classes = evalin('base','classes');
    [~, W_lda] = lda(X,classes);
    Xm = bsxfun(@minus, X, mean(X));
    Xproj = Xm * W_lda(:,1:2);
    set(handles.axes6,'visible', 'on')
    axes(handles.axes6);
    gscatter(Xproj(:,1), Xproj(:,2),classes);
    xL = xlim;
    yL = ylim;
    line([0 0],yL,'color', [0 0 0]);  %x-axis
    line(xL, [0 0],'color', [0 0 0]);  %y-axis
    title('LDA (original data)')
end

function elast_Callback(~, ~, ~)

function elast_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DTWButton_Callback(hObject, ~, handles)
asgn_fnames = handles.asgn_fnames;
brd = handles.brd;
elastval = get(handles.elast,'string');
elastval = floor(str2double(elastval));
handles.elastval = elastval;
guidata(hObject,handles);

if  size(handles.fnames(asgn_fnames(brd==1)),2)==2
    dtwdata = handles.reddata(:,handles.asgn_fnames(handles.brd == 1));
    one = dtwdata(:,1)';
    two = dtwdata(:,2)';
    [ds, ix, iy] = dtw(one,two,elastval);
    onewarp = one(:,ix);
    twowarp = two(:,iy);
    axes(handles.axes4)
    cla
    p4 = plot([1:1:length(onewarp)],onewarp,[1:1:length(twowarp)],twowarp);
    
    ylabel('amplitude')
    xlabel('profile depth')
    legend(handles.fnames(handles.asgn_fnames(handles.brd==1)));
elseif size(handles.fnames(asgn_fnames(brd==1)),2)==1 && get(handles.checkbox1,'value')==1 && isfield(handles, 'pdist');
    evaluatedNr = asgn_fnames(brd==1);
    edists = handles.pdist(:,evaluatedNr);
    [sedists, idx] = sort(edists);
    y= zeros(1,length(handles.fnames));
    y(:) = 1;
    y(idx(2:7))=2;
    y(idx(1))=3;
    
    axes(handles.axes4);
    cla reset
    p4 = gscatter(handles.Xproj(:,str2double(handles.pm1)),handles.Xproj(:,str2double(handles.pm2)),y, 'kgr');
    for i = 1:max(y)
        set(p4(i), 'Marker','s','MarkerSize', 6,'LineWidth', 1);
    end
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'color', [0 0 0]);  %x-axis
    line(xL, [0 0],'color', [0 0 0]);  %y-axis
    box on
    legend('DTW unalike','DTW alike','DTW template')
    handles.p4 = p4;
    guidata(hObject,handles)
    
else
    msgbox('Select a maximum of two datasets to compare DTW or enable DTW for displaying matches!')
end

function PCAButton_Callback(hObject, eventdata, handles)
axes(handles.axes4)
cla
plot(handles.redsteps(:,handles.asgn_fnames(handles.brd == 1)),handles.reddata(:,handles.asgn_fnames(handles.brd == 1)));
ylabel('amplitude')
xlabel('profile depth')
box on
set(handles.axes2,'visible','on')
legend(handles.fnames(handles.asgn_fnames(handles.brd==1)));
handles.gca = gca;
guidata(hObject,handles)



function cutXVal_Callback(~, ~, ~)
function cutXVal_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upLim_Callback(~, ~, ~)
function upLim_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in alignData.
function alignData_Callback(hObject, ~, handles)
reddata = handles.reddata;
redsteps = handles.redsteps;
steps = redsteps(:,1);
fnames = handles.fnames;
[~, stepZero]= min(abs(steps(steps<0)));
amp = '0';
while str2double(amp) < 1 || str2double(amp) > 60
    amp = inputdlg('Select amplitude [1-60%] where depth is 0','Data Alignment',1,{'30'});
end
amp = str2double(amp);
i = 0;
[startVal,shiftVal] = deal(zeros(1,length(fnames)));
for f = fnames
    i=i+1;
    startVal(i)= find(reddata(:,i) >= amp/100*max(reddata(:,i)),1);
    shiftVal(i) = startVal(i)-stepZero;
    reddata(:,i) = circshift(reddata(:,i),-shiftVal(i));
    if shiftVal(i)<0
        reddata(end-abs(shiftVal(i)):end,i)  =0;
    elseif shiftVal > 0
        reddata(1:shiftVal(i),i) = 0;
    end
    
end
axes(handles.view_data);
cla reset
colormap(dlmread('wb.map'));
imagesc(redsteps(:,1),1:1:size(reddata,2),reddata')
set(handles.view_data,'visible','on');
guidata(hObject,handles)


