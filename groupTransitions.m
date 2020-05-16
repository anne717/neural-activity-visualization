%Photometry analysis for SCN VIP sleep states- collating data over multiple mice: AV 18th April 2018
stage = 'REM';
stage2 = 'WAKE';
Acq_rate = 1; % user input (Hz): choose time resolution of final figure
pre_state_change_time = 60; % user input (seconds): choose time to display (and average) before a state change
post_state_change_time = 60; % user input (seconds): choose time to display after a state change
timeValuePre = 10; % user input (seconds): time want to analyse pre-transition for quantification
timeValuePost = 10;% user input (seconds): time want to analyse post-transition for quantification
timeZero = 0;
timeOff = 1200;
pre_state_time = pre_state_change_time*Acq_rate; %no longer in seconds, now in data points (so 600 data points for 1 minute if Acq_rate =10Hz)
post_state_time = post_state_change_time*Acq_rate; % in data points
timeValuePre = timeValuePre*Acq_rate; % in data points
timeValuePost = timeValuePost*Acq_rate;% in data points
cmin = -21; %value for yaxis (minimum) heatmap
cmax = 15; %value for yaxis (minimum) heatmap
ymin = -3.5;
ymax = 2.5;

%Select folder
folder_name = strcat(pwd,'\results\');

cd(folder_name)

names = dir(['vip*morning*_',stage,'to', stage2,'*.mat']);
filename =string;
groupedTrials = struct;

for i = 1:length(names)
    filename(i) = names(i).name;
end 

numMice = length(filename);

for i = numMice:-1:1    % Iterate backwards to fill the structure!
    field = char(filename(i));
    groupedTrials(i).mouseNames = field(1:end-4);
    groupedTrials(i).shortNames = field(1:6);
    trialStruct = load(filename(i));
    groupedTrials(i).trials = trialStruct.inverse_f_trials;
end

trialNums = zeros(numMice,1);

for i = 1:length(trialNums)
    groupedTrials(i).trials = groupedTrials(i).trials';
end

mean_groupedTrials = zeros(numMice,size(groupedTrials(i).trials,2));

for i = 1:length(trialNums)
    trialNums(i) = size(groupedTrials(i).trials,1);
    mean_groupedTrials(i,:) = mean(groupedTrials(i).trials,1);
end

cumTrialNums = cumsum(trialNums);
totalTrials = sum(trialNums);
lengthTrial = size(groupedTrials(1).trials,2);
trialsMatrix = zeros(totalTrials,lengthTrial);

startIndx = 0;

 for i = 1:length(trialNums)
    singleMouseTrials = groupedTrials(i).trials;
    stopIndx = size(groupedTrials(i).trials,1);
    trialsMatrix(startIndx+1:startIndx+stopIndx,:) = singleMouseTrials;
    startIndx = startIndx+stopIndx;
 end
 
 startIndx =0;
 minTrialNums = min(trialNums);
 trialsMatrixMin = zeros(minTrialNums*numMice,lengthTrial);
 
 for i = 1:length(trialNums)
    singleMouseTrials = groupedTrials(i).trials;
    stopIndx = minTrialNums;
    trialsMatrixMin(startIndx+1:startIndx+stopIndx,:) = singleMouseTrials(1:minTrialNums,:);
    startIndx = startIndx+stopIndx;
 end
 
meanMinPre = mean(trialsMatrixMin(:,pre_state_time-timeValuePre:pre_state_time),2);
meanMinPost = mean(trialsMatrixMin(:,pre_state_time:timeValuePost+pre_state_time),2);
meanMinPreMatrix = reshape(meanMinPre,[minTrialNums,numMice]);
meanMinPostMatrix = reshape(meanMinPost,[minTrialNums,numMice]);
 
flipTrialsMatrix = flipud(trialsMatrix); 
x = [(timeZero-pre_state_time) (timeZero+post_state_time)];
y = [totalTrials 1];
x_stateLine = [timeZero,timeZero];
y_stateLine = [1,totalTrials];
xOff_stateLine = [timeOff,timeOff];

shortNames={groupedTrials([1:end]).shortNames};
shortNames = string(shortNames);

condition =groupedTrials(2).mouseNames;
charCondition = char(condition);
condCondition = charCondition(8:end);

x_axis = [x(1):1:x(2)];
x_axis2 = linspace(0-pre_state_time,0+post_state_time, (pre_state_time+post_state_time)+1);
trialsMatrix2 = trialsMatrix(:,1:length(x_axis2));
flipTrialsMatrix2 = flipud(trialsMatrix2); 

figure(1)
hold();
imagesc(x,y,flipTrialsMatrix2)
plot(x_stateLine,y_stateLine,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
%plot(xOff_stateLine,y_stateLine,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
colormap parula
yticks(cumTrialNums)
yticklabels({shortNames})
h = colorbar;
ylabel(h,'\deltaF/F (%)')
set(get(h,'ylabel'),'rotation',-90)
set(get(h,'ylabel'),'position',[2.20 2.7280 0])
title(['Heatmap of GCaMP6s signal at ' condCondition ' in SCN^{VIP} neurons']);
ylabel('Mouse')
xlabel('time around a state change (seconds)')
hold();
saveas(gcf,strcat(folder_name, condCondition,'_all trials heatmap.png'));

y_stateLine2 = [min(trialsMatrix(:)) max(trialsMatrix(:))];
x_axis = [x(1):1:x(2)];
x_axis2 = linspace(0-pre_state_time,0+post_state_time, (pre_state_time+post_state_time)+1);
trialsMatrix2 = trialsMatrix(:,1:length(x_axis2));

figure(2)
hold;
plot(x_axis2, trialsMatrix2');
plot(x_stateLine,y_stateLine2,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
%plot(xOff_stateLine,y_stateLine2,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
title(['GCaMP6s signal at ' condCondition ' in SCN^{VIP} neurons']);
ylabel('\deltaF/F')
xlabel('time around a state change (seconds)')
hold;
saveas(gcf,strcat(folder_name, condCondition,'_all trials.png'));

figure(3)
hold;
[y_stateLine3] = plotStdErr(x_axis2,trialsMatrix2);
plot(x_stateLine,y_stateLine3,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
%plot(xOff_stateLine,y_stateLine3,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
title(['Average GCaMP6s signal at ' condCondition ' in SCN^{VIP} neurons over all mice']);
ylabel('\deltaF/F')
xlabel('time around a state change (seconds)')

% Set x-tick labels (so x-axis is in seconds)
xl = xticklabels;
xl = str2double(xl);
xl = xl/Acq_rate;
xticklabels({xl});
hold;
saveas(gcf,strcat(folder_name, condCondition,'_all trials avg.png'));

figure(4)
hold;
plot(x_axis2, mean_groupedTrials');
plot(x_stateLine,y_stateLine2,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
%plot(xOff_stateLine,y_stateLine2,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
title(['Mean GCaMP6s signal at ' condCondition ' in SCN^{VIP} neurons for each mouse']);
ylabel('\deltaF/F')
xlabel('time around a state change (seconds)')
hold;
saveas(gcf,strcat(folder_name, condCondition,'_all mouse avgs.png'));

figure(5)
hold;
[y_stateLine4] = plotStdErr(x_axis2,mean_groupedTrials);
plot(x_stateLine,y_stateLine4,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
%plot(xOff_stateLine,y_stateLine3,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
title(['Average of average GCaMP6s signals at ' condCondition ' in SCN^{VIP} neurons over all mice']);
ylabel('\deltaF/F')
xlabel('time around a state change (seconds)')

% Set x-tick labels (so x-axis is in seconds)
xl = xticklabels;
xl = str2double(xl);
xl = xl/Acq_rate;
xticklabels({xl});
hold;
saveas(gcf,strcat(folder_name, condCondition,'_avg of avgs.png'));

y_stateLine1 = [1,minTrialNums*numMice];
minTicks = (minTrialNums:minTrialNums:minTrialNums*numMice);
y1 = [1,minTrialNums*numMice];

x_axis = [x(1):1:x(2)];
x_axis2 = linspace(0-pre_state_time,0+post_state_time, (pre_state_time+post_state_time)+1);
trialsMatrixMin2 = trialsMatrixMin(:,1:length(x_axis2));
flipTrialsMatrixMin2 = flipud(trialsMatrixMin2); 

maxTrialsMatrixMin = max(flipTrialsMatrixMin2,[],2);
[E,index] = sortrows(maxTrialsMatrixMin);
[sortedFlipTrialsMatrixMin2] = flipTrialsMatrixMin2(index,:);
%sortedFlipTrialsMatrix2 = flipud(sortedFlipTrialsMatrix2);

clims = [cmin cmax];
figure(6)
hold();
%imagesc(x,y1,sortedFlipTrialsMatrixMin2)
imagesc(x,y1,sortedFlipTrialsMatrixMin2, clims)
plot(x_stateLine,y_stateLine1,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
%plot(xOff_stateLine,y_stateLine,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
colormap parula
yticks(minTicks)
yticklabels({shortNames})
h = colorbar;
ylabel(h,'\deltaF/F (%)')
set(get(h,'ylabel'),'rotation',-90)
set(get(h,'ylabel'),'position',[2.20 2.7280 0])
title(['Heatmap of GCaMP6s signal at ' condCondition ' in SCN^{VIP} neurons']);
ylabel('Mouse')
xlabel('time around a state change (seconds)')
hold();
saveas(gcf,strcat(folder_name, condCondition,'_MIN trials heatmap.png'));

figure(7)
hold;
[y_stateLine4] = plotStdErr(x_axis2,flipTrialsMatrixMin2);
plot(x_stateLine,y_stateLine4,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
ylim([ymin ymax])
%plot(xOff_stateLine,y_stateLine3,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
title(['Average GCaMP6 signal at ' condCondition ' in SCN^{VIP} neurons, equal trials/mouse']);
ylabel('\deltaF/F')
xlabel('time around a state change (seconds)')

% Set x-tick labels (so x-axis is in seconds)
xl = xticklabels;
xl = str2double(xl);
xl = xl/Acq_rate;
xticklabels({xl});
hold;
saveas(gcf,strcat(folder_name, condCondition,'_MIN trials avg.png'));