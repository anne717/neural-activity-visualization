function segmentFunction(mouseName, settings, export_folder);
%% Use this script to analyse set times pre and post transition state (WITHOUT UV
% correction)

%Photometry analysis for SCN VIP sleep states: AV 23rd July 2018
% Uses a 2nd order exponential fit to baseline both gcamp and uv signals
% DOES NOT use this to subtract background signal
% This is the deltaF/F used in the matrix??
% Includes only data before the light pulse (code 7)... though still need
% to 'not include' transitions in which post_time goes over the light pulse
%clear all

%%%% 16th May 2018: AV modified so that: 
% 1. Only transitions to state with a length as long as stageDur are analysed
% 2. Transitions including the LED noise are removed 
% 3. Transitions too close to the beginning or end of the recording are removed
% 4. Transitions from one stage (stage2) to another (stage) are analysed
% 5. Data includes the pre_state_time and the post_state_time (unlike
% sleepAV_functions3 and 4)
% 6. Transitions including the light pulse (7) are removed 
% 7. Uses 'baseline' to calculate gcamp_baselined- this CAN calculates the median
% valus for 60 secs pre-transition and subtracts this value from the whole transition
% event but I have disabled this so just taking the matrix around the
% transition points
%% 

%set expt variables
%folder_name = 'D:\VIP_photometry\november reexported files_noUV\sleepWake\longDur\';
%mouseName = 'vip360_evening'; 
%mouseName = names(i).name(1:end-4);
Mouse = strcat(mouseName,'.mat');

stage = settings.stage;     
stage2 = settings.stage2;                
stateDur = settings.stateDur; %length of time after the transition mouse must remain in state to include in analysis
statePreDur = settings.statePreDur; %length of time before the transition mouse must be in state to include in analysis
lightTag = settings.lightTag;
normTag = settings.normTag; %(0 for normal %, 1 for z-score, 2 for percentage of min and max point)
nanTag = settings.nanTag;  %if 1,  fitting gcamp_ds and uv_ds with NaNs)
Acq_rate = settings.Acq_rate; % user input (Hz): choose time resolution of final figure
pre_state_change_time = settings.pre_state_change_time; % user input (seconds): choose time to display (and average) before a state change
post_state_change_time = settings.post_state_change_time; % user input (seconds): choose time to display after a state change
timeZero = settings.timeZero; % enter time at which want vertical line drawn (usually at point of state change which is 0 seconds)

filter_tag = 0; % if 1, execute 'fix bleaching' commands

if stage == 1
    stageText = 'WAKE';
else if stage == 2
        stageText = 'NREM';
    else if stage == 3
            stageText = 'REM';
        else if stage == 8
                stageText = 'REM transition';
            end
        end
    end
end

if stage2 == 1
    stageText2 = 'WAKE';
else if stage2 == 2
        stageText2 = 'NREM';
    else if stage2 == 3
            stageText2 = 'REM';
        else if stage2 == 8
                stageText2 = 'REM transition';
            end
        end
    end
end

% Calculate timings
load(Mouse)
sampling_rate = 1/Calcium.interval; % Time intervals 

ds_factor = sampling_rate/Acq_rate; % calculate the downsampling factor
pre_state_time = pre_state_change_time*Acq_rate; %no longer in seconds, now in data points (so 600 data points for 1 minute if Acq_rate =10Hz)
post_state_time = post_state_change_time*Acq_rate; % in data points
time_window = pre_state_time + post_state_time+1; %in data points


%% Extract Calcium Signal (and do UV subtraction)
%get Ca2+ signal from loaded data
load(Mouse)

gcamp_signal = Calcium.values;
uv_signal = Backgrou.values;

state_times1 = instant_state.times;
codes1 = instant_state.codes(:,1);

%Extracting time data
N = Calcium.length; % Number of datapoints
dt = Calcium.interval; % Time intervals 
tmin = Calcium.start; % First Time value (in seconds)
tmax = Calcium.start + (N*dt); % Time at the end of data (in seconds)
t = linspace(tmin, tmax, N); % Equi-distant time vector (for FFT)

if lightTag ==1
    lightPulse = find(codes1 == 7); % find index of light pulse
    lightTime = state_times1(lightPulse);
    lightOnIndx = find(t>lightTime,1)-1;
    gcamp_signal = gcamp_signal(1:lightOnIndx);
    uv_signal = uv_signal(1:lightOnIndx);
    t = t(1:lightOnIndx);
else
end

%gcamp_signal = signalNorm;
gcamp_ds = group_z_project_vector(gcamp_signal, ds_factor);    %down sample to 'Acq_rate' by taking the mean every 10 samples (in our case)
uv_ds = group_z_project_vector(uv_signal, ds_factor);
time_ds = group_z_project_vector(t', ds_factor);

if nanTag == 1
    gcampValid = ~isnan(gcamp_ds);
    uvValid = ~isnan(uv_ds);
    f=fit(time_ds(gcampValid)', gcamp_ds(gcampValid)','exp2');   % exp2: Y = a*exp(b*x)+c*exp(d*x)
    f2=fit(time_ds(uvValid)', uv_ds(uvValid)','exp2');
else
    f=fit(time_ds', gcamp_ds','exp2');   % exp2: Y = a*exp(b*x)+c*exp(d*x)
    f2=fit(time_ds', uv_ds','exp2');
end

% Re-plotting the data...
%f=fit(time_ds', gcamp_ds','exp2');   % exp2: Y = a*exp(b*x)+c*exp(d*x)
%f2=fit(time_ds', uv_ds','exp2');    % poly: f(x) = p1*x + p2
                                    % exp1: f(x) = a*exp(b*x)

figure(1);
hold on
plot(f,'r',time_ds',gcamp_ds','b');
plot(f2,'r',time_ds',uv_ds','m');
hold off

coeffs1 = coeffvalues(f);
coeffs2 = coeffvalues(f2);

%fitData = coeffs1(1)*time_ds + coeffs1(2);
%fitData2 = coeffs2(1)*time_ds + coeffs2(2);

fitData = coeffs1(1)*exp(coeffs1(2)*time_ds) + coeffs1(3)*exp(coeffs1(4)*time_ds);
fitData2 = coeffs2(1)*exp(coeffs2(2)*time_ds) + coeffs2(3)*exp(coeffs2(4)*time_ds);

%fitData = coeffs1(1)*exp(coeffs1(2)*time_ds);
%fitData2 = coeffs2(1)*exp(coeffs2(2)*time_ds);

%figure(2)
%hold on
%plot(f,'r')
%plot(fitData,'g')
%hold off

normData = (gcamp_ds - fitData)./fitData;
normData2 = (uv_ds - fitData2)./fitData2;
%normData = (gcamp_ds(gcampValid) - fitData)./fitData;
%normData2 = (uv_ds(uvValid) - fitData2)./fitData2;
signalNorm2 = normData;
%signalNorm2 = normData-normData2; %uncomment this line to get background channel subtraction
minSignalNorm = min(signalNorm2);
maxSignalNorm = max(signalNorm2);
diffSignalNorm = maxSignalNorm - minSignalNorm;
stdData = nanstd(signalNorm2);
meanData = nanmean(signalNorm2);
zData = (signalNorm2-meanData)/stdData;

if normTag == 0
    signalNorm = signalNorm2*100;
else if normTag == 1
        signalNorm = zData;
    else if normTag == 2
            signalNorm =(signalNorm2-minSignalNorm)/diffSignalNorm*100;
        end
    end
end

%normData = (gcamp_ds - fitData)./fitData;
%normData = (gcamp_ds - fitData)+ median(gcamp_ds);
%normData2 = (uv_ds - fitData2)+ median(gcamp_ds);
%normData2 = (uv_ds - fitData2)./fitData2;
%signalNorm = normData-normData2; %uncomment this line to get background channel subtraction
%signalNorm = normData;

figure(3)
hold on
plot(time_ds, normData, 'b')
plot(time_ds, normData2, 'm')
plot(time_ds, signalNorm, 'r')
hold off

%fix bleaching
if filter_tag == 1
gcamp_ds_smooth = smooth(gcamp_ds,1000);
gcamp_fix = (gcamp_ds - gcamp_ds_smooth')+ median(gcamp_ds);
%gcamp_ds = gcamp_fix;
else
end

gcamp_norm = signalNorm;
%% Find code of the light pulse and don't go past this

%%get state change times for plots
codes = instant_state.codes(:,1);

lightPulse = find(codes == 7); % find index of light pulse
wake_indx = find(codes == stage); %find index of wake epochs

for i = 1:length(wake_indx)
    if wake_indx(i) > lightPulse
        wake_indx = wake_indx(1:i-1);    
    break
    end
end

%% Ensure we won't try and analyse data that goes beyond the beginning and ending of the recording
    
stateCount = 0;

if wake_indx(1)-1 <  1
    wake_indx = wake_indx(2:end);
end

if wake_indx(end)+1>length(codes)
    wake_indx = wake_indx(1:end-1);
end
%% 

%wake_indx = wake_indx(1:lightPulse-1); %only take state transitions before the light pulse 

%% Find what state mouse is in- if in correct state, include that transition in analysis
for i = 1:length(wake_indx)
    if codes(wake_indx(i)-1)== stage2
        stateCount = stateCount+1;
    end
end

selectedWake = zeros(stateCount,1);
j = 0;

for i = 1:length(wake_indx)
    if codes(wake_indx(i)-1)== stage2
        j = j+1;
       selectedWake(j)=(wake_indx(i));
    end
end
%% 

%% Find times of transitions in which the state is of at least length 'stateDur' and time spent in previous state is of at least 'statePreDur'
wake = nonzeros(selectedWake);
statePlusOneIndx = selectedWake +1;

%wake = wake_indx + 1; %uncomment this line to get transitions FROM a state
%for i = 1:length(wake)
 %   if wake(i) >  wake_indx(end)
  %  wake = wake(1:i-1);
   % end
%end
%wake = wake(2:end-1);
%sws = find(codes == 2);

state_times = instant_state.times; % gives times in seconds
wake_times = state_times(selectedWake);
wake_timesTrunc = state_times(selectedWake); %wake_times5 =  wake times that have pre and post state times that don't go past the beginning or ending of the recording 
statePlusOneTimes = state_times(statePlusOneIndx);
stateIndxTrunc = selectedWake;

k=0;

for i = 1:length(statePlusOneTimes)
    k = k+1;
    if statePlusOneTimes(i) - wake_times(i) < stateDur
      wake_timesTrunc = vertcat(wake_timesTrunc(1:(k-1)),wake_timesTrunc((k+1):end));
      stateIndxTrunc = vertcat(stateIndxTrunc(1:(k-1)),stateIndxTrunc((k+1):end));
      k = k-1;  
    end
end

stateMinusOneIndx = stateIndxTrunc - 1;
stateMinusOneTimes = state_times(stateMinusOneIndx);
state_timesTrunc = wake_timesTrunc;

k = 0;
for i = 1:length(stateMinusOneTimes)
    k = k+1;
    if wake_timesTrunc(i) - stateMinusOneTimes(i) < statePreDur
      state_timesTrunc = vertcat(state_timesTrunc(1:(k-1)),state_timesTrunc((k+1):end));
      stateIndxTrunc = vertcat(stateIndxTrunc(1:(k-1)),stateIndxTrunc((k+1):end));
      k = k-1;  
    end
end



wake_timesTrunc2 = state_timesTrunc;

for i = 1:length(wake_timesTrunc)
    if wake_timesTrunc(i)-pre_state_change_time < 1
        wake_timesTrunc2 = wake_timesTrunc2((i+1):end);
    else
        if wake_timesTrunc(i) + post_state_change_time > tmax
            wake_timesTrunc2 = wake_timesTrunc2(1:i-1);
            break
        end
    end
end
%% Find where this is noise in the recording and don't analyse transitions that have the noise in them (i.e. where noise falls in the pre or post_state_change_time)

noise_indx1 = find(codes == 9); %noise in the signal generated by the stimulation LED going off)
noise_indx2 = find(codes ==7); 
noise_indx = [noise_indx1, noise_indx2];
noise_indx = sort(noise_indx);
noise_times = state_times(noise_indx);

wake_timesFloor = floor(wake_timesTrunc2); % makes an integer by taking the floor
wake_timesIntegers = nonzeros(wake_timesFloor); % removes zero values 
wake_timesNoNoise = wake_timesIntegers;

k=0;

for i =1:length(wake_timesIntegers)
    k=k+1;
    for j = 1:length(noise_times)
    if  wake_timesIntegers(i)-(pre_state_change_time+2) < noise_times(j) && noise_times(j) < wake_timesIntegers(i)+(post_state_change_time+2)
        wake_timesNoNoise = vertcat(wake_timesNoNoise(1:(k-1)),wake_timesNoNoise((k+1):end));
        k = k-1;
    end
    end
end

stateChangeIndx = zeros(length(wake_timesNoNoise),1);

for i = 1:length(wake_timesNoNoise)
    stateChangeIndx(i) = find(time_ds>wake_timesNoNoise(i),1);
end
%% 
    
%wake_times3 = wake_times4;    
    
%wake_times3 = zeros(length(wake_times1),1);

%for i = 1:length(wake_times1)
 %   wake_times3(i) = find(time_ds>wake_times1(i),1);
%end

%wake_times1 = floor(wake_times2.*Acq_rate);


%sws_times = floor(state_times(sws));
A = zeros(length(gcamp_norm),1);
A(stateChangeIndx) =median(gcamp_norm);
A(A==0) = NaN;

figure(4);
%plot(gcamp_ds)
hold on
plot(gcamp_norm,'r')
hold on
% plot(photodi_ds*3.5,'g')
hold on
plot(A,'k-*')
title(['GCaMP6s signal with transition times'])

[gcamp_baselined,f0] = baseline(stateChangeIndx, time_window, pre_state_time, post_state_time, signalNorm);
uv_baselined = baseline(stateChangeIndx, time_window, pre_state_time, post_state_time, normData2);

inverse_uv = uv_baselined';
mean_uv = mean(inverse_uv,2);


x_axis = linspace(0-pre_state_time,0+post_state_time, (pre_state_time+post_state_time)+1);
figure(5);
inverse_f_trials = gcamp_baselined';
plot(x_axis, inverse_f_trials);
title(['GCaMP6s signal at transitions from ', stageText2, ' to ', stageText,' in SCN^{VIP} neurons: ', mouseName]);
saveas(gcf,strcat(export_folder,mouseName,'_', stageText2, 'to ', stageText,'_allTrials_noUV.png'));

% Calculating mean +/- SEM
mean_f_trials = mean(inverse_f_trials,2);
std_f_trials = std(inverse_f_trials,0,2);
div_factor = size(inverse_f_trials,2);
sem_f_trials = std_f_trials./sqrt(div_factor);
sem_plus = mean_f_trials + sem_f_trials;
sem_plus = sem_plus'; % need to invert 1 dimensional array so can peform calculations on it
sem_minus = mean_f_trials - sem_f_trials;
sem_minus = sem_minus'; % need to invert 1 dimensional array so can peform calculations on it

flip_sem_minus = fliplr(sem_minus); % need to flip sem_minus left to right to can draw the SEM polygon
sem_Y = [sem_plus,flip_sem_minus]; % get Y co-ordinates (in order) for SEM polygon
flip_x_axis = fliplr(x_axis); % flip x-axis left to right so can draw the SEM polygon
sem_x_axis = [x_axis,flip_x_axis]; % get X co-ordinates (in order) for SEM polygon


x_stateLine = [timeZero,timeZero];
y_stateLine = [min(sem_minus),max(sem_plus)];

figure(6);
hold on
plot(x_axis,mean_f_trials, 'r', 'LineWidth',2);
%plot(x_axis,mean_uv, 'b', 'LineWidth',2);
plot(x_axis,sem_plus,'r');
plot(x_axis,sem_minus, 'r'); 
plot(x_stateLine,y_stateLine,'Color',[0 0 1], 'LineWidth',1);
fill(sem_x_axis,sem_Y,'r','FaceAlpha',0.2, 'EdgeColor', 'none');
title(['GCaMP6s signal at transitions from ', stageText2, ' to ', stageText,' in SCN^{VIP} neurons: ', mouseName]);
ylabel('\deltaF/F')
xlabel('time around a state change (seconds)')

% Set x-tick labels (so x-axis is in seconds)
xl = xticklabels;
xl = str2double(xl);
xl = xl/Acq_rate;
xticklabels({xl});

hold off
saveas(gcf,strcat(export_folder,mouseName,'_', stageText2, 'to ', stageText,'noUV.png'));

%normTrials = inverse_f_trials- inverse_uv;
%normAvg = mean(normTrials,2);
x = [(timeZero-pre_state_time) (timeZero+post_state_time)];
y = [size(stateChangeIndx)];

figure(7)
hold();
imagesc(x,y,gcamp_baselined)
plot(x_stateLine,y,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
%plot(xOff_stateLine,y_stateLine,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
colormap parula
%yticks(cumTrialNums)
%yticklabels({shortNames})
h = colorbar;
ylabel(h,'\deltaF/F')
set(get(h,'ylabel'),'rotation',0)
title(['Heatmap of GCaMP6s signal at transitions from ', stageText2, ' to ', stageText,' in SCN^{VIP} neurons: ', mouseName]);
ylabel('Transitions')
xlabel('time around a state change (seconds)')
hold();
saveas(gcf,strcat(export_folder,mouseName,'_', stageText2, 'to ', stageText,'noUV_heatmap.png'));

%figure(8)
%plot(x_axis,normAvg, 'g', 'LineWidth',2);

save(strcat(export_folder,mouseName,'_', stageText2, 'to ', stageText,'noUV.mat'),'inverse_f_trials');
close all
end