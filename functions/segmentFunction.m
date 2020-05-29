function segmentFunction(mouseName, settings, export_folder)

%% Define settings
Mouse = strcat(mouseName,'.mat');
%Mouse = 'vip384_morning.mat';

stage = settings.stage;     % 1: WAKE, 2: NREM, 3: REM, 5: DOUBT
stage2 = settings.stage2;                
stateDur = settings.stateDur; %length of time after the transition mouse must remain in state to include in analysis
statePreDur = settings.statePreDur; %length of time before the transition mouse must be in state to include in analysis
normTag = settings.normTag; %(0 for normal %, 1 for z-score, 2 for percentage of min and max point)
plotFit = settings.plotFit;  %if 1,  save fitted curve
uvCorrectionTag = settings.uvCorrectionTag; % If 1, peform background subtraction 
Acq_rate = settings.Acq_rate; % user input (Hz): choose time resolution of final figure
pre_state_change_time = settings.pre_state_change_time; % user input (seconds): choose time to display (and average) before a state change
post_state_change_time = settings.post_state_change_time; % user input (seconds): choose time to display after a state change
timeZero = settings.timeZero; % enter time at which want vertical line drawn (usually at point of state change which is 0 seconds)

% Match stage codes in instant_state to sleep states
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

%% Load data from mouse
load(Mouse)

% Calculate time windows
pre_state_time = pre_state_change_time*Acq_rate; %no longer in seconds, now in data points (so 600 data points for 1 minute if Acq_rate =10Hz)
post_state_time = post_state_change_time*Acq_rate; % in data points
time_window = pre_state_time + post_state_time+1; %in data points

%get Ca2+ and UV signal from loaded data
gcamp_signal = Calcium.values;
uv_signal = Backgrou.values;

% Get sleep stage (state) marker codes and times
state_times = instant_state.times;
codes = instant_state.codes(:,1);

%Extracting timing data from Calcium struct
N = Calcium.length; % Number of datapoints
dt = Calcium.interval; % (Sampling interval from each data point in Calcium struct) 
tmin = Calcium.start; % First Time value (in seconds)
tmax = Calcium.start + (N*dt); % Time at the end of data (in seconds)
t = linspace(tmin, tmax, N); % Equi-distant time vector (for FFT)

% Determine if there is a light stimulus and only analyse up to this point
% if so

lightPulse = find(codes == 7); % find index of light pulse
if size(lightPulse)>0
  lightTag = 1;
else
    lightTag = 0;
end

if lightTag ==1
    lightTime = state_times(lightPulse);
    lightOnIndx = find(t>lightTime,1)-1;
    gcamp_signal = gcamp_signal(1:lightOnIndx);
    uv_signal = uv_signal(1:lightOnIndx);
    t = t(1:lightOnIndx);
else
end

% Determine if there are NaNs in the data
nanCheck = sum(isnan(Calcium.values));
if nanCheck>0
    nanTag = 1;
else
    nanTag = 0;
end

%down sample Calcium and UV signal to 'Acq_rate' by taking the mean every 'ds_factor' samples (in our case)
ds_factor = (1/dt)/Acq_rate;
gcamp_ds = group_z_project_vector(gcamp_signal, ds_factor);    
uv_ds = group_z_project_vector(uv_signal, ds_factor);
time_ds = group_z_project_vector(t', ds_factor);

%% Fit the Calcium and UV data to 2nd order exponential curves
if nanTag == 1
    gcampValid = ~isnan(gcamp_ds);
    uvValid = ~isnan(uv_ds);
    f=fit(time_ds(gcampValid)', gcamp_ds(gcampValid)','exp2');   % exp2: Y = a*exp(b*x)+c*exp(d*x)
    f2=fit(time_ds(uvValid)', uv_ds(uvValid)','exp2');
else
    f=fit(time_ds', gcamp_ds','exp2');   % exp2: Y = a*exp(b*x)+c*exp(d*x)
    f2=fit(time_ds', uv_ds','exp2');
end

% Plot Calcium and UV signal, together with fitted curves
figure(1);
hold on;
plot(f,'r',time_ds',gcamp_ds','b');
plot(f2,'r',time_ds',uv_ds','m');
hold off;

if plotFit == 1
    saveas(gcf,strcat(export_folder,mouseName,'_fittedData.png'));
end

% Calculate coefficients of the fitted curves
coeffs1 = coeffvalues(f);  % Calcium coefficients
coeffs2 = coeffvalues(f2);  %UV coefficients

% Re-create curves using the same timebase as Calcium and UV signals s cna
% perform dF/F calculation
fitData = coeffs1(1)*exp(coeffs1(2)*time_ds) + coeffs1(3)*exp(coeffs1(4)*time_ds);
fitData2 = coeffs2(1)*exp(coeffs2(2)*time_ds) + coeffs2(3)*exp(coeffs2(4)*time_ds);

% Calculate dF/F
if nanTag == 1
    normData = (gcamp_ds(gcampValid) - fitData(gcampValid))./fitData(gcampValid);
    normData2 = (uv_ds(uvValid) - fitData2(uvValid))./fitData2(uvValid);
else
    normData = (gcamp_ds - fitData)./fitData;
    normData2 = (uv_ds - fitData2)./fitData2;
end

if uvCorrectionTag == 0
    signalNorm2 = normData;
else 
    signalNorm2 = normData-normData2;
end

%% Normalize data based on normTag
minSignalNorm = min(signalNorm2);
maxSignalNorm = max(signalNorm2);
diffSignalNorm = maxSignalNorm - minSignalNorm;
stdData = nanstd(signalNorm2);
meanData = nanmean(signalNorm2);
zData = (signalNorm2-meanData)/stdData;

if normTag == 0
    signalNorm = signalNorm2*100;
else
    if normTag == 1
        signalNorm = zData;
    else
        if normTag == 2
            signalNorm =(signalNorm2-minSignalNorm)/diffSignalNorm*100;
        end
    end
end

% Plot dF/F for Calcium and UV, together with normalized data
figure(3)
hold on;
if nanTag == 1
    plot(time_ds (gcampValid), normData, 'b')
    plot(time_ds(uvValid), normData2, 'm')
    signalNormValid = ~isnan(signalNorm);
    plot(time_ds(signalNormValid), signalNorm, 'r')     
else
    plot(time_ds, normData, 'b')
    plot(time_ds, normData2, 'm')
    plot(time_ds, signalNorm, 'r')
end 
hold off;


%% Get state change times for plots
gcamp_norm = signalNorm;
wake_indx = find(codes == stage); %find index of wake epochs

% If there is a light stimulus, only analyse data up until this point
for i = 1:length(wake_indx)
    if lightTag ==1
        if wake_indx(i) > lightPulse
            wake_indx = wake_indx(1:i-1);    
        break
        end
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

%% Find what state mouse is in- if in correct state, include that transition in analysis
stateCount = 0;
for i = 1:length(wake_indx)
    if codes(wake_indx(i)-1)== stage2
        stateCount = stateCount+1;
    end
end

if stateCount<1
    message = ['No transitions from ' stageText2  ' to ' stageText ' in ' mouseName];
    disp(message)
    close all
    return
else
end

    
selectedWake = zeros(stateCount,1);
j = 0;

for i = 1:length(wake_indx)
    if codes(wake_indx(i)-1)== stage2
        j = j+1;
       selectedWake(j)=(wake_indx(i));
    end
end

%% Find times of transitions in which the state is of at least length 'stateDur' and time spent in previous state is of at least 'statePreDur'
wake = nonzeros(selectedWake);
statePlusOneIndx = selectedWake +1;

wake_times = state_times(selectedWake);
wake_timesTrunc = state_times(selectedWake); 
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

for i = 1:length(state_timesTrunc)
    if state_timesTrunc(i)-pre_state_change_time < 1
        wake_timesTrunc2 = wake_timesTrunc2((i+1):end);
    end
end

wake_timesTrunc2 = state_timesTrunc;
for i = 1:length(state_timesTrunc)
    if state_timesTrunc(i) + post_state_change_time > tmax
        wake_timesTrunc2 = wake_timesTrunc2(1:i-1);
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

%% Plot normalized data with asterisks showing where state transitions occur
A = zeros(length(gcamp_norm),1);
A(stateChangeIndx) =median(gcamp_norm);
A(A==0) = NaN;

figure(4);
hold on;
plot(gcamp_norm,'r')
plot(A,'k-*')
title(['GCaMP6s signal with transition times']);
hold off;

%% Segment data into a matrix where each row contains the pre/post time data for a transiton
[gcamp_baselined,f0] = baseline(stateChangeIndx, time_window, pre_state_time, post_state_time, signalNorm);
uv_baselined = baseline(stateChangeIndx, time_window, pre_state_time, post_state_time, normData2);

% Construct x-axis
x_axis = linspace(0-pre_state_time,0+post_state_time, (pre_state_time+post_state_time)+1);

% Plot data for each transition overlaid
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

% Plot mean+/ SEM data over all transitions
figure(6);
hold on;
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

hold off;
saveas(gcf,strcat(export_folder,mouseName,'_', stageText2, 'to ', stageText,'noUV.png'));


x = [(timeZero-pre_state_time) (timeZero+post_state_time)];
y = [size(stateChangeIndx)];

% Plot each transition as a heatmap
figure(7)
hold on;
imagesc(x,y,gcamp_baselined)
plot(x_stateLine,y,'Color',[0 0 0], 'LineWidth',1, 'LineStyle', '--');
colormap parula
h = colorbar;
ylabel(h,'\deltaF/F')
set(get(h,'ylabel'),'rotation',0)
title(['Heatmap of GCaMP6s signal at transitions from ', stageText2, ' to ', stageText,' in SCN^{VIP} neurons: ', mouseName]);
ylabel('Transitions')
xlabel('time around a state change (seconds)')
hold off;
saveas(gcf,strcat(export_folder,mouseName,'_', stageText2, 'to ', stageText,'noUV_heatmap.png'));

% Save matrix containing segmented transition data into export data for
% future analysis
save(strcat(export_folder,mouseName,'_', stageText2, 'to ', stageText,'noUV.mat'),'inverse_f_trials');
close all

end