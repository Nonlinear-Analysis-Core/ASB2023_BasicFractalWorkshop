%% Basic Fractal Analysis in Movement Science Workshop - ASB 2023 

% DESCRIPTION:
% This code is intended to get you familiar with creating your own time
% series for simulations It also guides you through the logic behind some
% of the parameter selections for DFA as well as being a guide on how to
% use the NONAN DFA function in your own research. We provide a finally 
% fully worked example on how to implement the function in your analysis
% pipeline.
    % Section 1 - Pink, white and brown noise
    % Section 2 - DFA and time series length
    % Section 3 - Changing DFA Parameters
    % Section 4 - DFA analysis on a single walking timeseries
    % Section 5 - DFA analysis on a single COP timeseries
    % Section 6 - DFA analysis on multiple walking timeseries
    % Section 7 - DFA analysis on multiple COP timeseries


%% Section 1: Pink, white and brown noise

% As we have learnt, the alpha value of a time series can vary depending on
% what the time series actually is. In the section we will create some time
% series with a pink, white and brown power spectrum. This code recreates 
% one of the figures in the workshop presentation.

% Parameters for DFA
order = 1;
n_min = 10;
n_max = floor(2000/8);
n_length = 18;
plotOption = 1;

figure()
% Anti-correlated Time Series
tiledlayout(2,4)
nexttile
hurst_exp = 0.2;
anti_corr = pinkNoise(hurst_exp, 1500);
plot(anti_corr, 'k', 'LineWidth', 1)
xlabel('Observations')
ylabel('Amplitude (a.u.)')
title('Anti-Correlated Noise')

% White Noise
nexttile
hurst_exp = 0.5;
white = pinkNoise(hurst_exp, 1500);
plot(white,'k', 'LineWidth', 1)
xlabel('Observations')
ylabel('Amplitude (a.u.)')
title('White Noise')

% Pink Noise
nexttile
hurst_exp = 0.99;
pink = pinkNoise(hurst_exp, 1500);
plot(pink,'k', 'LineWidth', 1)
xlabel('Observations')
ylabel('Amplitude (a.u.)')
title('Pink Noise')

% Brown Noise
nexttile
hurst_exp = 0.5;
brown = cumsum(pinkNoise(hurst_exp, 1500));
plot(brown,'k', 'LineWidth', 1)
xlabel('Observations')
ylabel('Amplitude (a.u.)')
title('Brownian Motion')

% Anti-correlated Time Series DFA plot
nexttile
[a, r2, out_a, out_l] = dfa(anti_corr, n_min, n_max, n_length, plotOption);

% White noise DFA plot
nexttile
[a, r2, out_a, out_l] = dfa(white, n_min, n_max, n_length, plotOption);

% Pink noise DFA plot
nexttile
[a, r2, out_a, out_l] = dfa(pink, n_min, n_max, n_length, plotOption);

% Brown noise DFA plot
nexttile
[a, r2, out_a, out_l] = dfa(brown, n_min, n_max, n_length, plotOption);

set(gcf, 'Position', [50 80 1300 650])

%%
% Now we have created some times series to play around with and to look at
% other statistical properties of. Lets use pink and white noise as our two
% examples. 

% Pink Noise
figure()
tiledlayout(2,2)
nexttile([1 2])
plot(pink, 'k', 'LineWidth', 1)

nexttile
histogram(pink, 'FaceColor', 'k', 'FaceAlpha', 0.7)

nexttile
autocorr(pink, 50)
xlim([0 50])

% White Noise
figure()
tiledlayout(2,2)
nexttile([1 2])
plot(white, 'k', 'LineWidth', 1)

nexttile
histogram(white, 'FaceColor', 'k', 'FaceAlpha', 0.7)

nexttile
autocorr(white, 50)
xlim([0 50])

%% Section 2 - DFA and time series length

% In this section we will load in the summary results of some simulations.
% To show that the length of the time series matters, we performed DFA on
% time series' of various lengths and expected alpha values. 2000
% simulations were performed for each combination of length and expected
% alpha value.

% The graphs show us a number of interesting things. Firstly, the DFA is
% positively biased - especially when the time series is really short.
% Panel two shows as the time series length increased, the absolute
% difference between the actual outcome and the expected alpha decreases
% dramatically until we hit 500 data points and then it starts to plateau.
% The same trend can be seen in the standard deviation of the calculated 
% alpha values. This demonstrates the importance of having enough data 
% points in your time series to reliably estimate alpha.

T = readtable("SimStats.csv"); % Load in the summary results from the simulations

% Plot the mean of the calculated alpha vs the expected alpha value
tiledlayout(3,1);
nexttile
scatter(T.Length, T.Actual_Mean, [], T.Expected, 'filled')
yline(T.Expected, '--k', 'LineWidth', 0.5)
ylabel('Mean of Alpha')
title('Alpha Mean vs Time series length')
set(gca, 'FontSize', 15)
xlim([180 1120])
xticks([])
ylim([0 1.1])

% Plot the absolute difference between calculated alpha vs the expected alpha value
nexttile
scatter(T.Length, T.diff, [], T.Expected, 'filled')
ylabel('Abs. Difference of Alpha')
title('Abs. Difference of Actual Alpha values vs Expected Values')
xlim([180 1120])
xticks([])
set(gca, 'FontSize', 15)

% Plot the SD of the calculated alpha vs the expected alpha value
nexttile
scatter(T.Length, T.Actual_SD, [], T.Expected, 'filled')
xlabel('Time Series Length')
ylabel('SD of Alpha')
title('SD vs Time series length')
xlim([180 1120])
set(gca, 'FontSize', 15)

% Setting the colourbar properties and also the size of the plot
cb = colorbar;
cb.Ticks = unique(T.Expected);
cb.TickLabels = [0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
cb.Layout.Tile = 'east';
cb.Label.String = 'Alpha Value';
set(gca, 'FontSize', 15)
set(gcf, 'Position', [250 150 1000 700])

%% Section 3: Changing DFA Parameters

% During the workshop we talked about recommendations for some paramters...
% let's see what happens to the resulting alpha values when we change these
% parameters a little. Commented out below is the code to replicate the
% simulations that we performed. You can run section if you would like but
% it does take quite a while to finish (>2hrs). We recommend changing the
% sims value to something smaller if you are going to run on a laptop.

%% Simulations for changing n_min
% sims = 1000;
% ts_length = 1000;
% hurst = 0.1:0.1:1;
% hurst(10) = 0.99;
% 
% % Changing n_min
% n_min = 8:2:27;
% n_max = ts_length/8; % We will use a length of 1000 so it is hard coded here. 
% n_length = 18; % default number of points to sample best fit.
% plotOption = 1; % return a plot
% 
% res = [];
% 
% updates = 1:20:1000; % Loop numbers that you want to save on
% 
% for s = 1:sims
% 
%     for l = 1:length(n_min)
%         disp(n_min(l))
% 
%         for h = 1:length(hurst)
% 
%             this = pinkNoise(hurst(h), ts_length);
%             a(l,h) = dfa(this, n_min(h), n_max, n_length, plotOption);
%             H(h) = hurst(h);
% 
%         end % Hurst loop
% 
%         res = [res; repmat(n_min(l), 10,1), a(l,:)', H'];
% 
%     end % Time series length loop
% 
%     if ismember(s, updates)
%         writematrix(res, 'n_minSims2.csv');
%         disp("Saving data.....")
%     else
%         % Do nothing
%         fprintf("Continuing on with loop %d.\n", s)
% 
%     end % End if statement
% 
% end % Sim number loop
% 
% T = array2table(res);
% T.Properties.VariableNames = {'n_min','Actual','Expected'};
% writetable(T, 'n_minSims2.csv');
% 
% 
% nminSims = readtable("n_minSims.csv");
% 
% vals = unique(nminSims.n_min);
% final = [];
% for i = 1:length(unique(nminSims.n_min))
% 
%     % I want to take the length, expected and actual columns from the group
%     % stats and add them to the bottom of an array. 
% 
%     j = nminSims(nminSims.n_min == vals(i),:);
%     j1 = grpstats(j, ["n_min", "Expected"], ["mean", "std"]);
% 
%     final = [final; j1.n_min, j1.Expected, j1.mean_Actual, j1.std_Actual];
% 
% end
% T = array2table(final);
% %T.Properties.VariableNames = {'Length','0.1','0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', 'Hurst'};
% T.Properties.VariableNames = {'n_min','Expected', 'Actual_Mean', 'Actual_SD'};
% difference = abs(T.Actual_Mean - T.Expected);
% T.diff = difference;
% 
% writetable(T, "n_minStats.csv")

%% Changing DFA Parameters Continued.....

% replicate for regression
% Exercise: Randomly delete some data and see how DFA reacts 

n_minStats = readtable("n_minStats.csv");

% Plot the n_min used vs the mean alpha value
tiledlayout(3,1)
nexttile
scatter(n_minStats.n_min, n_minStats.Actual_Mean, [], n_minStats.Expected, 'filled')
yline(n_minStats.Expected, '--k', 'LineWidth', 0.5)
ylabel('Mean of Alpha')
title('n_min vs Alpha Mean')
xlim([7.5 26.5])
xticks([])
set(gca, 'FontSize', 15)

% Plot the n_min used vs the absolute difference between calculated alpha
nexttile
scatter(n_minStats.n_min, n_minStats.diff, [], n_minStats.Expected, 'filled')
ylabel('Abs. Difference of Alpha')
title('Abs. Difference of Actual Alpha values vs Expected Values')
xlim([7.5 26.5])
xticks([])
set(gca, 'FontSize', 15)

% Plot the n_min used vs the SD of the actual alpha values
nexttile
scatter(n_minStats.n_min, n_minStats.Actual_SD, [], n_minStats.Expected, 'filled')
xlabel('# of data points in smallest window')
ylabel('SD of Alpha')
title('SD vs Time series length')
xlim([7.5 26.5])
set(gca, 'FontSize', 15)

% Setting the colourbar properties and also the size of the plot
cb = colorbar;
cb.Ticks = unique(n_minStats.Expected);
cb.TickLabels = [0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
cb.Layout.Tile = 'east';
cb.Label.String = 'Alpha Value';
set(gca, 'FontSize', 15)
set(gcf, 'Position', [250 150 1000 700])


% EXERCISE: 
%   - Can you write some code to change the n_length parameter? Does
%     changing the number of points in the regression line change the alpha
%     value returned?

%% Section 4: DFA analysis on one stride interval timeseries

% Now that we have gotten a good grasp of the basics with simulated time
% series, we will jump in with some real data! 
% This section of the tutorial will briefly show how to run DFA on stride
% interval data from a single participant while walking at their self-selected
% pace overground. 

% Clearing everything to start from scratch
clear all;
close all;
clc;

% Loading in data. Column 1 is time, Column 2 is stride intervals.
ts = readmatrix('S308_SPW_ex2.csv');

% ---------- Defining our input parameters for the dfa.m function ---------

% n_min is set to 16 as recommended. The citation can be found on line 25
% of the dfa.m function.
n_min = 16;

% n_max set to a default of the length of the time series divided by 9
% as recommended by Damouras, et al., 2010. The appropriate paper can be
% found on line 25 of the dfa.m function.
n_max = length(ts)/9;

% n_length is set to the default number of points to sample best fit.
n_length = 18;

% Sets the plotOption command to true so we return a plot.
plotOption = 1;

% The next line of code runs the actual DFA analysis. This function runs
% through the steps outlined in the lecture. The output, a, represents the
% alpha value and is what we are most interested in when comparing between
% subjects, trials, conditions etc. Take note of how the alpha value
% compares to the DFA analysis interpretation. Is it persistent or
% antipersistent?

% Panel 1: Stride Intervals overtime
tiledlayout(2,3)
nexttile([1 3])
plot(ts, 'k', 'LineWidth', 1.0);
title('DFA Stride Intervals: S208');
xlabel('Time (s)');
ylabel('Stride Interval (seconds)');
set(gca, 'FontSize', 14); % Setting Font Size

% Panel 2: DFA (Exact same process as Section 1)
nexttile
[a, r2] = dfa(ts, n_min, n_max, n_length, plotOption);
title('DFA Stride Intervals');
xlabel('log(n)');
ylabel('log(F(n))');
set(gca, 'FontSize', 14); % Setting Font Size

% Panel 3: Histogram of stride intervals overtime
nexttile
histogram(ts, 'FaceColor', 'k', 'FaceAlpha', 0.7);
title('Stride Interval Distribution');
xlabel('Stride Interval (seconds)');
ylabel('Frequency of Occurance');
set(gca, 'FontSize', 14); % Setting Font Size

% Panel 4: Autocorrelation of stride intervals 
nexttile
autocorr(ts, 50);
title('Stride Interval Autocorrelation');
xlim([0 50])
xticks(0:10:50)
set(gca, 'FontSize', 14); % Setting Font Size
set(gcf, 'Position', [300, 100, 1100, 800]); % Setting Figure Size


%% Section 5: DFA analysis on one COP timeseries

% Center of Pressure data for Subject 28 of the Human Balance Evaluation
% Database (https://archive.physionet.org/physiobank/database/hbedb/)
% The subject performed a 60 second standing task with their eyes closed
% on a firm surface.

% This section of the tutorial will briefly show how to run DFA on center
% of pressure velocity for a single subject for a single trial. Notice,
% this procedure is the exact same as Section 1.

% Clearing everything to start from scratch
clear all;
close all;
clc;

% Loading in data. Column 1 is time, Column 3 is COPy Velocity.
s018_cop = readmatrix('s018_cop.csv');

% ----------------- Defining our input parameters for DFA -----------------

% Select column 3 from our data which is COPy velocity
ts = s018_cop(:,3);

% n_min is set to 16 as recommended
n_min = 16;

% n_max set to a default of the length of the time series divided by 9
% as recommended by Damouras, et al., 2010.
n_max = length(ts)/9;

% n_length is set to the default number of points to sample best fit.
n_length = 18;

% Sets the plotOption command to true so we return a plot.
plotOption = 1;

% The next line of code runs the actual DFA analysis. This function runs
% through the steps outlined in the lecture. The output, a, represents the
% alpha value and is what we are most interested in when comparing between
% subjects, trials, conditions etc. Is it persistent or antipersistent?

[a, r2, out_a, out_l] = dfa(ts, n_min ,n_max, n_length, plotOption);
title('DFA COPy Velocity s018');
xlabel('log(n)');
ylabel('log(F(n))');
set(gca, 'FontSize', 18); % Setting Font Size
set(gcf, 'Position',  [800, 100, 1000, 800]); % Setting Figure Size;


% The regression line doesn't really capture the trend of the data. To get
% around this we will find where the data curves and perform a regression
% either side of that breakpoint. 
breakpoint = 98;

n = out_a(:,1); % box sizes
F = out_a(:,2); % fluctuation for a given box size
logn = log(n); % log box sizes
logF = log(F); % log fluctuation for a given box size

fontSize = 15;
markerSize = 20;
[ashort,along] = crossover(logn, logF, breakpoint, 1);
title('DFA COPy Velocity short and long fluctuation');
xlabel('log(n)');
ylabel('log(F(n))');
caption = sprintf(' long \\alpha = %.4f\n short \\alpha = %.4f ', along, ashort);
text(min(logn), max(logF), caption, 'FontSize', fontSize);

% As we can see, dividing the data at the break point gives us a much
% better fit of the line. What does it mean though? Can you think of a
% way to automate this process?

%% Section 6: DFA analysis on multiple stride interval timeseries

% This data was collected as part of an experiment where participants
% were asked to synchronise their heel strikes to different pacing signals.
% Raffalt PC et al. DOI: 10.1016/j.neulet.2022.136909

% Clearing everything to start from scratch
clear all;
close all;
clc;

% Setup for the loop

% Gets directory to run the rest of the code. Simply choose where
% the data is stored. In this case it is in the same directory as the
% script.
my_directory = uigetdir;

% We want to store the results in a specific results folder so we
% will create that here.
output_directory = append(my_directory, '/ANALYSIS OUTPUT');

% This is a little unnecessary but shows that you can create a folder from
% within MATLAB. Something like this may be useful if you want to create
% many folders in a loop and want to name them using a subject ID.
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% Because all our files are in the same directory we can use the wildcard
% '*' character and part of a string that we want to find to choose what
% data we will use.
my_files = dir(fullfile(my_directory, '*ex2.csv'));

% Create empty array to store trial information in
output = [];


for k = 1:length(my_files)
    base_filename = my_files(k).name;
    full_filename = fullfile(my_directory, base_filename);

    % Loading in data
    ts = readmatrix(full_filename);

    % Get the file name. This will be used to create the figure title and
    % also the .png file names.
    filename_split = strsplit(base_filename, '_');

    % ------ Defining our input parameters for the dfa.m function -------

    % n_min is set to 16 as recommended
    n_min = 16;

    % n_max set to a default of the length of the time series divided by 9
    % as recommended by Damouras, et al., 2010.
    n_max = length(ts)/9;

    % n_length is set to the default number of points to sample best fit.
    % This input can be moved outside the loop if preferred.
    n_length = 18;

    % Sets the plotOption command to true so we return a plot.
    plotOption = 1;

    % Up until this point the code is not much different than Section 1. Now we
    % will produce a single figure with our time series plotted, our DFA
    % output, and a histogram of our time series.

    % --------------- Run DFA and create a plot of results ---------------

    % Figure 1: Stride Intervals overtime
    tiledlayout(2,2)
    nexttile([1 2])
    plot(ts, 'k', 'LineWidth', 1.0);
    title(filename_split{1},' ', filename_split{2}]);
    xlabel('Time (s)');
    ylabel('Stride Interval (seconds)');
    set(gca, 'FontSize', 18); % Setting Font Size

    % Figure 2: DFA (Exact same process as Section 1)
    nexttile
    [a, r2] = dfa(ts, n_min, n_max, n_length, plotOption)
    title('DFA Stride Intervals');
    xlabel('log(n)');
    ylabel('log(F(n))');
    set(gca, 'FontSize', 18); % Setting Font Size

    % Figure 3: Histogram of stride intervals overtime
    nexttile
    histogram(ts);
    title('Stride Interval Distribution');
    xlabel('Stride Interval (seconds)');
    ylabel('Frequency of Occurance');
    set(gca, 'FontSize', 18); % Setting Font Size
    set(gcf, 'Position', [100, 100, 1000, 800]); % Setting Figure Size

    % The next two lines of code will require you to press any key on the
    % keyboard to close the figure and progress to the next trial
    disp('Press any key to continue!')
    pause;

    % -------------------- Data collating for analysis --------------------

    % Collating our figures into the folder titled "ANALYSIS OUTPUT"
    png = append(filename_split{1}, '.png');
    saveas(gcf,fullfile(output_directory, png));

    % Collating alpha values from our DFA function with corresponding
    % subject id and condition.
    output{k,1} = filename_split{1};
    output{k,2} = filename_split{2};;
    output{k,3} = a;

    close all % Close figures before looping through to the next trial
end

% Exporting DFA values into a table in long format that can be easily used
% for statistical analysis and even more visualisations.
dfa_results = cell2table(output, 'VariableNames', {'id', 'condition', 'alpha'});
writetable(dfa_results, fullfile(output_directory, 'Stride Interval DFA Values.csv'));

%% Section 7: DFA analysis on multiple COP timeseries

% Clearing everything to start from scratch
clear all;
close all;
clc;

% Setup for the loop

% Gets directory to run the rest of the code. Simply choose where
% the data is stored. In this case it is in the same directory as the
% script.
my_directory = uigetdir;

% We want to store the results in a specific results folder so we
% will create that here.
output_directory = append(my_directory, '\ANALYSIS OUTPUT');

% This is a little unnecessary but shows that you can create a folder from
% within MATLAB. Something like this may be useful if you want to create
% many folders in a loop and want to name them using a subject ID.
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% Because all our files are in the same directory we can use the wildcard
% '*' character and part of a string that we want to find to choose what
% data we will use.

my_files = dir(fullfile(my_directory, '*cop.csv'));

dfa_values = []; % This is where our DFA reults will go

for k = 1:length(my_files)
    base_filename = my_files(k).name;
    full_filename = fullfile(my_directory, base_filename);

    % Loading in data. Column 1 is time, Column 3 is COPy Velocity
    temp_dat = readmatrix(full_filename);

    % -------------- Defining our input parameteres for DFA --------------

    % Select column 2 from our data corresponding to velocity
    ts = temp_dat(:,3);

    % Get the file name. This will be used to create the figure title and
    % also the .png file names.
    filename_split = strsplit(base_filename, '.');

    % n_min is set to 16 as recommended.
    n_min = 16;

    % n_max set to a default of the length of the time series divided by 9
    % as recommended by Damouras, et al., 2010.
    n_max = length(ts)/9;

    % n_length is set to the default number of points to sample best fit.
    % This input can be moved outside the loop if preferred.
    n_length = 18;

    % Sets the plotOption command to true so we return a plot.
    plotOption = 1;

    % Up until this point the code is not much different than Section 1. Now we
    % will produce a single figure with our time series plotted, our DFA
    % output, and a histogram of our time series.

    % ------------------------- Plot and run DFA -------------------------

    % Figure 1: COPy velocity overtime
    subplot(2, 2, [1,2]);
    plot(ts);
    title(filename_split{1}, 'Interpreter', 'none');
    xlabel('Time (s)');
    ylabel('Velocity (cm/s)');
    set(gca, 'FontSize', 18); % Setting Font Size

    % Figure 1: DFA (Exact same process as Section 1)
    subplot(2, 2, 3);
    [a, r2] = dfa(ts, n_min, n_max, n_length, plotOption)
    title('DFA COPy Velocity');
    xlabel('log(n)');
    ylabel('log(F(n))');
    set(gca, 'FontSize', 18); % Setting Font Size

    % Figure 3: Histogram of COPy velocity overtime
    subplot(2, 2, 4);
    histogram(ts);
    title('Velocity Distribution');
    xlabel('Velocity (cm/s)');
    ylabel('Frequency of Occurance');
    set(gca, 'FontSize', 18); % Setting Font Size
    set(gcf, 'Position', [800, 100, 1000, 800]); % Setting Figure Size;

    % -------------------- Data collating for analysis --------------------

    % Collating our figures into the folder titled "ANALYSIS OUTPUT"
    png = append(filename_split{1}, '.png');
    saveas(gcf,fullfile(output_directory, png));

    % Collating alpha values from our DFA function with corresponding
    % subject id and condition.
    dfa_values{k,1} = filename_split{1};
    dfa_values{k,2} = a;

    close all % Close figures before looping through to the next trial
end

% Exporting DFA values into a table in long format that can be easily used
% for statistical analysis and even more visualisations.
dfa_results = cell2table(dfa_values, 'VariableNames', {'id','alpha'});
writetable(dfa_results, fullfile(output_directory, 'COP DFA Values.csv'));

% Done!
