%% Lab 6 Neural Analysis II
% Student Name: Vanessa Ma
clear
warning('off');
load('Lab6_CenterOutTest.mat');
load('ML.mat');
load('PV.mat');

trials = length(go);
N = length(unit);

%%% Data Cleaning
W = zeros(N, trials);
time_buckets = sort([go - 1; go + 1]);
angle_buckets = 0: 45: 315;
spikes_arr = zeros(N, trials);
for n = 1:N
    if ~isempty(unit(n).times)
        trials_clean = histcounts(unit(n).times, time_buckets);
        spikes_arr(n, :) = trials_clean(1:2:trials*2);
    else
        spikes_arr(n, :) = 0;
    end
    % basic
%     W(n, :) = spikes_arr(n, :);
    % QUESTION 3: baselining method 1
    W(n, :) = spikes_arr(n, :) - mean(spikes_arr(n, :));
    % QUESTION 3: baselining method 2
%     W(n, :) = spikes_arr(n, :) - tuning_arr(n, 1);
end

%%% Population Vector Testing (EXERCISE 1)

PV_x = zeros(1, N);
PV_y = zeros(1, N);

for n = 1:N
    PV_x(n) = cos(pref_arr(n));
    PV_y(n) = sin(pref_arr(n));
end

test_x = zeros(1, trials);
test_y = zeros(1, trials);

for trial = 1:trials
    test_x(trial) = PV_x * W(:, trial);
    test_y(trial) = PV_y * W(:, trial);
end

dir_array = mod(180 * atan2(test_y, test_x) / pi, 360);
pv_ans = zeros(2, trials);
for dir = 1:8
    if dir == 1
        bucket_idx = find(dir_array > 315 + 22.5 & dir_array <= 22.5);
    else
        bucket_idx = find(dir_array > angle_buckets(dir) - 22.5 & dir_array <= angle_buckets(dir) + 22.5);
    end
    pv_ans(1, bucket_idx) = dir;
end

% accuracy of the method (EXERCISE 3 & QUESTION 1)
pv_ans(2, :) = boolean(pv_ans(1,:) == direction');
pv_p_correct = sum(pv_ans(2,:)) / trials;
fprintf('Accuracy of Population Vector Method: %d \n', pv_p_correct);

%%% Maximum Likelihood Testing (Poisson)
p_spike = zeros(N, dir);
p_log = zeros(trials, dir);
ml_ans= zeros(2, trials);
for trial = 1:trials
    for dir = 1:8
        % poisson (EXERCISE 2)
        p_spike(:, dir) = poisspdf(spikes_arr(:, trial), mean_r(:, dir));
        % gaussian (QUESTION 4)
%         p_spike(:, dir) = normpdf(spikes_arr(:, trial), mean_r(:, dir), sd_r(:, dir));
        p_idx = find(~isnan(p_spike(:, dir)));
        p_clean = p_spike(p_idx, dir);
        p_log(trial, dir) = sum(log(p_clean));
    end
    [ml_mx, ml_ans(1, trial)] = max(p_log(trial, :));
end

% accuracy of the method (EXERCISE 3 & QUESTION 1)
ml_ans(2, :) = boolean(ml_ans(1,:) == direction');
ml_p_correct = sum(ml_ans(2,:)) / trials;
fprintf('Accuracy of Population Vector Method: %d \n', ml_p_correct);

% QUESTION 2 - Distribution of Preferred directions
figure;
histogram(direction, 'BinEdges', 0.5:1:8.5);
title('Predicted Directions for 80 Test Trials (Actual)');
xlabel('Preferred Direction');
ylabel('Number of Trials');
    % expected uniform distribution
% e_dist = unifrnd(1, 8, N, 1);
% e = trials/8*ones(8, 1);
% [h, p, stats] = chi2gof(direction, 'Expected', e, 'Edges', 0.5:1:8.5);
% fprintf('p-value of actual directions: %d \n', p);

figure;
histogram(ml_ans(1,:), 'BinEdges', 0.5:1:8.5);
title('Predicted Directions for 80 Test Trials (Maximum Likelihood, Poisson)');
xlabel('Preferred Direction');
ylabel('Number of Trials');
    % chi-squared test for maximum likelihood
% [h_ml, p_ml, stats_ml] = chi2gof(ml_ans(1,:), 'Expected', e, 'Edges', 0.5:1:8.5);
% fprintf('p-value of maximum likelihood predictions: %d \n', p_ml);
    
figure;
histogram(pv_ans(1,:), 'BinEdges', 0.5:1:8.5);
title('Predicted Directions for 80 Test Trials (Population Vector, Baselined)');
xlabel('Preferred Direction');
ylabel('Number of Trials');
    % chi-squared test for population vector
% [h_pv, p_pv, stats_pv] = chi2gof(pv_ans(1,:), 'Expected', e, 'Edges', 0.5:1:8.5);
% fprintf('p-value of population vector predictions: %d \n', p_pv);

% QUESTION 2 - Fit of cosine curve

