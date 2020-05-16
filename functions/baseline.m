function [All_f_trials,f0] = baseline1(state_time, window, pre, post, signal);
All_f_trials = zeros(length(state_time),window);
f0 = zeros(length(state_time));
for i = 1:length(state_time)
    if (((state_time(i)))-pre) < 0     %if state change happens less than 1 time window into the recording...
        for j = 1:window     
            All_f_trials(i,j) = 0;    %...make from -time_window to state change a 0
        end
    else if (((state_time(i)))+post)>size(signal,2)     %if state change c;happens less than 1 time window into the recording...
        for j = 1:window     
            All_f_trials(i,j) = 0;    %...make from -time_window to state change a 0   
        end
    else  
        f0(i) = mean(signal(state_time(i)-pre:state_time(i)-1)); %takes the mean of the values in the preceeding pre_state_time values to state change
        f1(i) = median(signal(state_time(i)-pre:state_time(i)-1));
        for k = 1:window
            %All_f_trials(i,k) = ((signal(state_time(i)-pre-1+k))-f0(i))/f0(i); % deltaF /F(t)  = F(t) - mean(F)/mean(F) 
            All_f_trials(i,:) = signal(state_time(i)-pre:state_time(i)+post)-f1(i); % deltaF /F(t)  = F(t) - mean(F)/mean(F) 
            %All_f_trials(i,:) = signal(state_time(i)-pre:state_time(i)+post);
        end
        end
    end
end