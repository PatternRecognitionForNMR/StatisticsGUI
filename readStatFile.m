%% read data real smooth
% readStatFile('C:\Users\Chris\ownCloud\MATLAB\Statistics\allprofiles-2.mat')
function [data_raw, fnames, alldata, steps] = readStatFile(str1,swNormalize)
path = str1; % nothing to see here, passing on the string normally fails in R2016b
[~,~,ext] = fileparts(path); % determine file extension

switch ext % file extension
    case '.mat'
        data_raw = load(path);
        i=0;
        fnames = fieldnames(data_raw)'; %mat file variables are now fieldnames
        [maxStep, minStep, stepSize, stepZero] = deal(zeros(1,length(fnames))); %preallocate
        %% find min and max for steps, also deal with stepsize
        for f = fnames
            i=i+1;
            maxStep(i) = max(data_raw.(f{1})(:,1));
            minStep(i) = min(data_raw.(f{1})(:,1));
            stepSize(i) = data_raw.(f{1})(2,1)-data_raw.(f{1})(1,1);
            if stepSize(i) < 0 % if negative
                data_raw.(f{1}) = flipud(data_raw.(f{1})); %stepsizes are inverted for better displaying later on
                stepSize(i) = stepSize(i)*-1;
                stepZero(i) = min(abs(data_raw.(f{1})(:,1)));
            end
        end
        %min and max for steps are found
        globMinStep = min(minStep);
        globMaxStep = max(maxStep);
        if all(stepSize == stepSize(1)) == 0 || all(stepZero == stepZero(1)) == 0  %check for uniform stepsize
            i =0;
            msgbox('Uh-Oh. Step Size is not uniform throughout all your Profiles. Interpolation is used, be sure to reduce the Number of Steps manually now.')
            for f = fnames % interpolate otherwise
                i=i+1;
                bSteps = minStep(i):1:maxStep(i);
                bData = interp1(data_raw.(f{1})(:,1),data_raw.(f{1})(:,2),minStep(i):1:maxStep(i));
                data_raw.(f{1}) = [bSteps ;bData]';
                clear stepSize
                stepSize = 1;
                
            end
        else
            stepSize(2:end) = [];
        end
        
        %% align the data
        steps = repmat(globMinStep:stepSize:globMaxStep,length(fnames),1)';
        alldata = zeros(size(steps));
        i =0;
        for f = fnames % interpolate otherwise
            i=i+1;
            bMember = ismember(steps(:,i), data_raw.(f{1})(:,1)); % find out where to sort in our data
            firstMention = find(bMember==1,1);
            lastMention = find(bMember==1,1,'last');
            k =0;
            for j = firstMention:1:lastMention
                k = k+1;
                alldata(j,i) =data_raw.(f{1})(k,2);
            end
        end
        
    case {'.txt','.dat','.csv','.xlsx'}
        tab = readtable(path); %use table nomenclature for simplicity, sacrificing speed, sorry winXP computers
        fnames = tab.Properties.VariableNames; % Profile names extracted from table along with row one name 'steps'
        fnames(1) = []; % delete variable name for steps
        
        data_raw = table2array(tab);
        steps = repmat(data_raw(:,1),1,length(fnames)); %this may be unnecessary
        alldata = data_raw(:,2:end);
        
        
    otherwise
        msgbox('Input file type not recognized. We support .dat, .csv, .xls, .txt and .mat');
end
%normalize data if wanted
if swNormalize ==1
    for i = 1:size(alldata,2)
        alldata(:,i)= alldata(:,i)/max(alldata(:,i));
    end
end

end