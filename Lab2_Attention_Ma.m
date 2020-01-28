%% LAB 2 - ATTENTION
%  Student Name: Vanessa Ma
%% Experimental Data Collection

% Global Variables
block = 4;
rounds = 4;
data = zeros(block*rounds, 8);
% set random number generator to be dependent on the cpu clock
rng('shuffle');
validity_arr = rand(block, 1);
validity = -1;
response = 0;
rt = 0;
response_corr = 0;

figure;
for j = 1:rounds
    for i = 1:block
        % fixation point
        plot(13, 13, 'k+');
        set(gca, 'XColor', 'none', 'YColor', 'none');
        xlim([0 26])
        ylim([0 26])
        right_far = rectangle('Position', [21,11,4,4]);
        left_far = rectangle('Position', [1,11,4,4]);
        right_close = rectangle('Position', [14,11,4,4]);
        left_close = rectangle('Position', [8,11,4,4]);
        box_array = [right_far, right_close, left_far, left_close];
        key_array = ['k', 'j', 'd', 'f'];
        ascii_array = [107, 106, 100, 102];
        instruct_key = zeros(1, length(key_array));
        % initialize variables
        instruct_text = 0;
        position = [];
        if(i == 1)
            for k = 1:length(key_array)
                position = get(box_array(k), 'Position');
                instruct_text = text(position(1) + 2, position(2) + 2, key_array(k));
                instruct_key(k) = instruct_text;
            end
            instruct_1 = text(6, 20, 'Press the following keys for the target square.');
            instruct_2 = text(8, 18, 'Press no key if no target appears.');
            instruct_3 = text(9, 9, 'Press any key to continue');
            pause;
            delete([instruct_key, instruct_1, instruct_2, instruct_3])
        end
        % ensure 15% neutral trials
        if(validity_arr(i) <= 0.15)
            validity = -1;
        % ensure 42.5% valid trials 
        elseif((validity_arr(i) > 0.15) && (validity_arr(i) <= 0.575))
            validity = 0;
        elseif(validity_arr(i) > 0.575)
            validity = 1;
        end
        pause(0.2);
        [response, rt, validity, response_corr, cue_coord, target_coord] = stim(validity, box_array, ascii_array);
        data((j-1)*block + i, :) = [response, rt, validity, response_corr, cue_coord, target_coord];
        pause(0.1);
    end
end
text(12, 14, 'End');

pause(0.3);
save('Lab2_data.mat', 'data')

%% Data Analysis

load('Lab2_data.mat')

% Drop all incorrect trials
data_clean = data(data(:,4) ==  1, :);

% Subset for all conditions
data_valid = data_clean(data_clean(:,3) == 1, :);
data_invalid = data_clean(data_clean(:,3) == 0, :);
data_neutral = data_clean(data_clean(:,3) == -1, :);

% Average reaction times and standard errors for valid and invalid condition
valid_mean = mean(data_valid(:, 2)); 
valid_se = std(data_valid(:, 2)) / sqrt(length(data_valid(:, 2)));
invalid_mean = mean(data_invalid(:, 2)); 
invalid_se = std(data_invalid(:, 2)) / sqrt(length(data_invalid(:, 2)));

% Create new column of cue and target distance
dist = zeros(length(data_invalid(:,1)), 1);
for i = 1:length(data_invalid(:,1))
    dist(i) = data_invalid(i, 5) - data_invalid(i, 7);
end
data_invalid = [data_invalid dist];

% Create array of means and corresponding reaction times
dist_array = [0, 6, 7, 13];
means_array = [valid_mean, 0, 0, 0];
for i = 2:length(means_array)
    dist_subset = data_invalid(data_invalid(:, 9) == dist_array(i), :);
    means_array(i) = mean(dist_subset(:, 2));
end

% Error Rates
%   valid trials
%   invalid trials
%   neutral trials

%%% WRITE-UP

%% Question 1: Plot the cue, target, fixation point on the same graph for an invalid trial
% copied code from stim.m
    % choose box for cue and valid target
v_side = randi([1,4]);
i_side = ceil(rand*4);
    % choose box for invalid target
while(v_side == i_side)
    i_side = ceil(rand*4);
end
figure;
plot(13, 13, 'k+');
set(gca, 'XColor', 'none', 'YColor', 'none');
xlim([0 26])
ylim([0 26])
right_far = rectangle('Position', [21,11,4,4]);
left_far = rectangle('Position', [1,11,4,4]);
right_close = rectangle('Position', [14,11,4,4]);
left_close = rectangle('Position', [8,11,4,4]);
box_array = [right_far, right_close, left_far, left_close];
box_array(v_side).LineWidth = 1;
box_array(v_side).EdgeColor = 'r';
target_pos = get(box_array(i_side), 'Position');
target = text(target_pos(1) + 2, target_pos(2) + 2,'*', 'FontSize', 24);

% Question 2

%%% Question 2a: Mean reaction time for valid versus invalid trials? 
disp('Mean reaction time, valid: ');
disp(valid_mean);
disp('Mean reaction time, invalid: ')
disp(invalid_mean);

%%% Question 2b: Is it statistically significant? 
% Significance Testing
[significant, p] = ttest2(data_valid(:, 2), data_invalid(:, 2));
disp('Statistical Significance (1/0 - Yes/No): ');
disp(significant);
disp('P-value: ');
disp(p);

%%% Question 2c: How do you expect varying the ratio of valid to invalid trials would affect the results?
% If the ratio of invalid trials to valid trials were increased, subjects
% would learn to expect that the target would not be in the cued box. 


% Question 3

%%% Question 3a: Plot the mean reaction time versus distance between target and cue (ranging from 0 to 3). 
figure;
% scatter(data_invalid(:, 2), data_invalid(:, 9));
scatter(dist_array, means_array);

%%% Question 3b: Estimate the slope of the line. 
line_x = linspace(min(dist_array), max(dist_array), 100);
means_fit = polyfit(dist_array, means_array, 1);
means_y = polyval(means_fit, line_x);
hold on
plot(line_x, means_y, 'b-');
hold off
[means_coef, means_p] = corrcoef(dist_array, means_array);

disp(means_coef);

%%% Question 3c: What is the speed of the 'attentional spotlight' in targets/second? What is the approximate speed in cm/s?
% The first question is a crude way of quantifying how many supposed
% 'targets' can be held in mind at the same time. The 
% The speed of the attentional spotlight in cm/s is simply the slope of the
% plot with distance on the y-axis and reaction time on the x-axis. (i.e.
% the reciprocal of the slope derived in 3b.
% 

% An interesting future experiment to run would be to validate the
% 'inhibition of return' conclusion in attention, where exogenous cues as
% that used in this experiment can briefly enhance (100-300ms ISI) speed
% and accuracy of a valid trial, but then impairs reaction time results for
% ISI beyond that point (500-3000ms).