names = dir('vip*_morning.mat');
settings.stage = 1;     % 1: WAKE, 2: NREM, 3: REM, 5: DOUBT
settings.stage2 = 3;     
           
settings.stateDur = 10; %length of time after the transition mouse must remain in state to include in analysis
settings.statePreDur = 10; %length of time before the transition mouse must be in state to include in analysis
settings.lightTag = 0;
settings.normTag = 1; %(0 for normal %, 1 for z-score, 2 for percentage of min and max point)
settings.nanTag = 1;  %if 1,  fitting gcamp_ds and uv_ds with NaNs)
settings.Acq_rate = 1; % user input (Hz): choose time resolution of final figure
settings.pre_state_change_time = 60; % user input (seconds): choose time to display (and average) before a state change
settings.post_state_change_time = 60; % user input (seconds): choose time to display after a state change
settings.timeZero = 0; % enter time at which want vertical line drawn (usually at point of state change which is 0 seconds)

current_folder = pwd;
if ~exist('results', 'dir')
    mkdir('results')
end

export_folder = strcat(current_folder,  '\results\');

for i = 1:length(names)
    mouseName = names(i).name(1:end-4);
    segmentFunction(mouseName, settings, export_folder);
end