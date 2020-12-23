%% LAB 1 VISUAL SEARCH
% Student Name: Vanessa Ma

%% Experimental Data Collection

% Global Variables Set Up
set_size_arr = [4 6 8 12];
block = 50;
% initialize data collection
data_p = zeros(block * length(set_size_arr), 4);
data_c = zeros(block * length(set_size_arr), 4);
% initialize target conditions
    % randomize trial presence
target_arr = randi([0,1], 1, block); 
    % target set up
target_col = 'r';
out_col = 'b';
target_char = 'x';
out_char = 'o';

% experiment procedure begins
for j = 1:length(set_size_arr)
    % Pop-out Block Instructions
    plot(50, 50, strcat(target_col, '+'));
    xlim([0 100])
    ylim([0 100])
    set(gca, 'XColor', 'none', 'YColor', 'none');
    text(20, 60, 'Press "1" if target color appears, "0" otherwise');
    text(30, 40, 'Press any key to continue');
    pause;
    for i = 1:block
        % Pop-out Trial
        target = target_arr(i);
        num = set_size_arr(j);
        clf;
        xlim([0 100])
        ylim([0 100])
        stim_plot_p_col(num, target, target_char, target_col, out_char, out_col);
        set(gca, 'XColor', 'none', 'YColor', 'none')
        % get response time
        tic;
        pause;
        time = toc;
        % get response
        response = get(gcf, 'CurrentCharacter');
        % check response correctness
        response_corr = 0;
        if(isequal(response, 48) && (target == 0))
            response_corr = 1;
        elseif(isequal(response, 49) && (target == 1))
            response_corr = 1;
        end
        % data_collection
        data_p((j-1)*block + i,:) = [num, response, time, response_corr];
        pause(0.2);
        clf;
    end
    text(0.3, 0.5, 'Press any key for next block');
    set(gca, 'XColor', 'none', 'YColor', 'none')
    pause;
    % Conjunction Block Instructions
    clf;
    plot(50, 50, strcat(target_col, target_char));
    xlim([0 100])
    ylim([0 100])
    set(gca, 'XColor', 'none', 'YColor', 'none');
    text(10, 60, 'Press "1" if target character and color appears, "0" otherwise');
    text(30, 40, 'Press any key to continue');
    pause;
    for i = 1:block
        % Conjunction Trial
        target = target_arr(i);
        num = set_size_arr(j);
        clf;
        xlim([0 100])
        ylim([0 100])
        stim_plot_c(num, target, target_char, target_col, out_char, out_col);
        set(gca, 'XColor', 'none', 'YColor', 'none')
        % get response time
        tic;
        pause;
        time = toc;
        % get response
        response = get(gcf, 'CurrentCharacter');
        % check response correctness
        response_corr = 0;
        if(isequal(response, 48) && (target == 0))
            response_corr = 1;
        elseif(isequal(response, 49) && (target == 1))
            response_corr = 1;
        end
        % data_collection
        data_c((j-1)*block + i,:) = [num, response, time, response_corr];
        pause(0.2);
    end
    clf;
    xlim([0 100])
    ylim([0 100])
    text(30, 50, 'Press any key for next block');
    set(gca, 'XColor', 'none', 'YColor', 'none')
    pause;
end
clf;
xlim([0 100])
ylim([0 100])
text(30, 50, 'End');
set(gca, 'XColor', 'none', 'YColor', 'none')
pause;
close;

% Saving data. 
save('Lab1_data.mat', 'data_p', 'data_c');

%% Data Analysis
% p - pop-out
% c - conjunction
% corr - correct trials only
% pos - target existent
% neg - target nonexistent

load('Lab1_data.mat');
% reload data variables for running this section and advance
set_size_arr = [4 6 8 12];

% Drop Incorrect Trials
data_p_corr = data_p(data_p(:,4) == 1, :);
data_c_corr = data_c(data_c(:,4) == 1, :);

% Understand error rates
%   Pop-out, by set size
p_four_accuracy = 1 - sum(sum(data_p_corr == 4)) / 50;
p_six_accuracy = 1 - sum(sum(data_p_corr == 6)) / 50;
p_eight_accuracy = 1 - sum(sum(data_p_corr == 8)) / 50;
p_twelve_accuracy = 1 - sum(sum(data_p_corr == 12)) / 50;
%       Print error rates
fprintf('Pop-Out Error Rates (out of 50 trials) \n')
fprintf('Set Size 4: %f \n', p_four_accuracy);
fprintf('Set Size 6: %f \n', p_six_accuracy);
fprintf('Set Size 8: %f \n', p_eight_accuracy);
fprintf('Set Size 12: %f \n', p_twelve_accuracy);
%   Conjunction, by set size
c_four_accuracy = 1 - sum(sum(data_c_corr == 4)) / 50;
c_six_accuracy = 1 - sum(sum(data_c_corr == 6)) / 50;
c_eight_accuracy = 1 -sum(sum(data_c_corr == 8)) / 50;
c_twelve_accuracy = 1 - sum(sum(data_c_corr == 12)) / 50;
%       Print error rates
fprintf('Conjunction Error Rates (out of 50 trials) \n')
fprintf('Set Size 4: %f \n', c_four_accuracy);
fprintf('Set Size 6: %f \n', c_six_accuracy);
fprintf('Set Size 8: %f \n', c_eight_accuracy);
fprintf('Set Size 12: %f \n', c_twelve_accuracy);

% Separate into positive and negative trials
data_p_corr_pos = data_p_corr(data_p_corr(:,2) == 49,:);
data_p_corr_neg = data_p_corr(data_p_corr(:,2) == 48,:);
data_c_corr_pos = data_c_corr(data_c_corr(:,2) == 49,:);
data_c_corr_neg = data_c_corr(data_c_corr(:,2) == 48,:);

% Mean and Standard Error for Reaction Times
    % initialize collection arrays
    % means - collect average reaction times; se - collect standard errors
p_pos_means = zeros(0,4);
p_pos_se = zeros(0,4);
p_neg_means = zeros(0,4);
p_neg_se = zeros(0,4);
c_pos_means = zeros(0,4);
c_pos_se = zeros(0,4);
c_neg_means = zeros(0,4);
c_neg_se = zeros(0,4);

for j = 1:length(set_size_arr)
    % For Pop-Out, With Targets
    p_pos_subset = data_p_corr_pos(data_p_corr_pos(:,1) == set_size_arr(j), :);
    p_pos_means(j) = mean(p_pos_subset(:,3));
    p_pos_se(j) = std(p_pos_subset(:,3)) / sqrt(length(p_pos_subset(:,3)));
    % For Pop-Out, Without Targets
    p_neg_subset = data_p_corr_neg(data_p_corr_neg(:,1) == set_size_arr(j), :);
    p_neg_means(j) = mean(p_neg_subset(:,3));
    p_neg_se(j) = std(p_neg_subset(:,3)) / sqrt(length(p_neg_subset(:,3)));
    % For Conjunction, With Targets
    c_pos_subset = data_c_corr_pos(data_c_corr_pos(:,1) == set_size_arr(j), :);
    c_pos_means(j) = mean(c_pos_subset(:,3));
    c_pos_se(j) = std(c_pos_subset(:,3)) / sqrt(length(c_pos_subset(:,3)));
    % For Conjunction, Without Targets
    c_neg_subset = data_c_corr_neg(data_c_corr_neg(:,1) == set_size_arr(j), :);
    c_neg_means(j) = mean(c_neg_subset(:,3));
    c_neg_se(j) = std(c_neg_subset(:,3)) / sqrt(length(c_neg_subset(:,3)));
end

% Linear Fit
line_x = linspace(min(set_size_arr), max(set_size_arr), 100);
    % For Pop-Out, With Targets
p_pos_fit = polyfit(set_size_arr, p_pos_means, 1);
p_pos_y = polyval(p_pos_fit, line_x);
    % For Pop-Out, No Targets
p_neg_fit = polyfit(set_size_arr, p_neg_means, 1);
p_neg_y = polyval(p_neg_fit, line_x);
    % For Conjunction, With Targets
c_pos_fit = polyfit(set_size_arr, c_pos_means, 1);
c_pos_y = polyval(c_pos_fit, line_x);
    % For Conjunction, No Targets
c_neg_fit = polyfit(set_size_arr, c_neg_means, 1);
c_neg_y = polyval(c_neg_fit, line_x);

%% Write-Up

% Global Variables Set Up for following code
target_col = 'r';
out_col = 'b';
target_char = 'x';
out_char = 'o';

%% Question 1: Experimental Procedure Plots
disp('1a) Pop-out Plot (set size = 12)');
figure;
xlim([0 100])
ylim([0 100])
stim_plot_p_col(12, 1, target_char, target_col, out_char, out_col);
disp('1b) Conjunction Plot (set size = 12)');
figure;
xlim([0 100])
ylim([0 100])
stim_plot_c(12, 1, target_char, target_col, out_char, out_col);

%% Question 2: Plot the 4 traces of mean reaction times vs. set size for all 4 trial types
% scatterplots of mean reaction times
figure;
scatter(set_size_arr, p_pos_means, 'b','s');
xlim([0 15])
xlabel('Set Size')
ylim([-0.1 0.9])
ylabel('Reaction Time (s)')
hold on
scatter(set_size_arr, p_neg_means, 'c', 's')
hold on
scatter(set_size_arr, c_pos_means, 'r', '*')
hold on
scatter(set_size_arr, c_neg_means, 'm', '*')
hold on
% lines of fit
plot(line_x, p_pos_y, 'b-')
hold on
plot(line_x, p_neg_y, 'c-')
hold on
plot(line_x, c_pos_y, 'r-')
hold on
plot(line_x, c_neg_y, 'm-')
legend('pop-out with target', 'pop-out without target', 'conjunction with target', 'conjunction without target', 'Location', 'northwest');
hold off

%% Question 3: Report correlation coefficient and p_values for each of the 4 trial types
% Pop-out, Target present
[p_pos_coef, p_pos_p] = corrcoef(set_size_arr, p_pos_means);
disp('Pop-Out, Target Present - Coefficients');
disp(p_pos_coef);
disp('Pop-Out, Target Present - P-Values');
disp(p_pos_p)
% Pop-out, Target not present
[p_neg_coef, p_neg_p] = corrcoef(set_size_arr, p_neg_means);
disp('Pop-Out, Target Not Present - Coefficents');
disp(p_neg_coef);
disp('Pop-Out, Target Not Present - P-Values');
disp(p_neg_p);
% Conjunction, Target present
[c_pos_coef, c_pos_p] = corrcoef(set_size_arr, c_pos_means);
disp('Conjunction, Target Present - Coefficients');
disp(c_pos_coef);
disp('Conjunction, Target Present - P-Values');
disp(c_pos_p);
% Conjunction, Target not present
[c_neg_coef, c_neg_p] = corrcoef(set_size_arr, c_neg_means);
disp('Conjunction, Target Not Present - Coefficients');
disp(c_neg_coef);
disp('Conjunction, Target Not Present - P-Values');
disp(c_neg_p);

%% Question 4

%%% Question 4a: How do the reaction times for pop-out compare to the 
%%% conjunction search when target is present?

% Comparing the blue and red lines, we can see that the reaction times for
% pop-out are generally lower than those for the conjunction search
% paradigm. As set size increases, the reaction time for pop-out conditions
% also remains roughly constant, while those for conjunction conditions
% scale somewhat linearly.

%%% Question 4b: For both searches, what effect does target presence have
%%% on reaction times?

% Target presence seems to lead to lower reaction times (i.e. faster
% response) in both cases.

%%% Question 4c: Is your data consistent with Treisman's conclusions?

% Roughly, yes. 1) Less time is needed to process singular features than
% multiple (conjunction) of features to separate items from their
% foreground, 2) singular feature processing is flat with set size, but 
% conjunction processing scales with stimuli complexity, 3) negative searches 
% generally take longer than target present searches. However, the lines of 
% best fit only loosely reflect the data and have some low correlation 
% coefficients (proxy for r-squared value for assessing fit), some down to 50%.
% None of the p-values were significant also to the 0.05 mark.

% Theres is also one deviation - the mean reaction time for lack of stimulus
% decreased as the set size increased for pop-out searches.

% However, insofar as Treisman's conclusion was that the human mind registers features of
% visual stimuli are perceived in a parallel fashion, and integrated in
% later stages of processing. From an experimental perspective, this means
% that feature identification (i.e. a pop-out search) should be far more
% efficient than conjunction identification, because the latter engages
% higher order synthesis functions. Broadly, these experimental results
% exhibit the same experimental conclusion.

% Some obvious limitations that may be the cause of the deviation are 1)
% low sample size and 2) possible bias, since I was both the experiment
% designer and the experimenter. Had I not known that the only two options
% in a block are pop-out target present / not present, I may have taken
% longer to parse the negative pop-out searches, as an example. Interesting
% further experiments to run would be character pop-out searches where
% stimuli could vary in color (or are all the same color) to understand if
% there are some more easily processed features. Alternatively, instead of
% specifying what the entire block is going to be, it might be interesting
% to run a mixed block with the same target (e.g. a red x). Lastly, noting
% error rates for conjunction tests are generally higher than those of
% pop-out tests, it may be worthwhile trying out a larger block size of
% conjunction trials in the future to ensure that average evaluation of
% both types of searches are truly comparable.

                
