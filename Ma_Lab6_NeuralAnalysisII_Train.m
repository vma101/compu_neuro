clear
warning('off');
load('Lab5_CenterOutTrain.mat');

trials = length(go);
N = length(unit);

mean_r = zeros(N, 8);
sd_r = zeros(N, 8);
fr = zeros(trials, N, 8);
pref_arr = zeros(1, N);

% fitting parameters
ang_fit = 0: 0.001: 2*pi;
cos_fun = @(p, theta) p(1) + p(2) *cos(theta - p(3));
tuning_arr = zeros(N, 3);
res_arr = zeros(N, 8);
fit_arr = zeros(N, length(ang_fit));


for n = 1:N
    % filter to only neurons that have spiking activity
    if ~isempty(unit(n).times)
        for dir = 1:8
            % find trials for each direction
            trial_idx = find(direction == dir);
            for trial = 1:length(trial_idx)
                cue = go(trial_idx(trial));
                % center spike times
                trial_rel = unit(n).times - cue;
                % extract only +/- 1s spikes
                trial_clean = trial_rel(trial_rel < 1 & trial_rel > -1);
                % compute firing rates for single trial
                fr(trial_idx(trial), n, dir) = length(trial_clean) / 2;
            end
            % data for MAXIMUM LIKELIHOOD method:
                % mean, sd of firing rates in all directions
            mean_r(n, dir) = mean(nonzeros(fr(:, n, dir)));
            sd_r(n, dir) = std(nonzeros(fr(:, n, dir)));
        end
    else
        mean_r(n, dir) = 0;
        sd_r(n, dir) = 0;
    end
    % fit training for POPULATION VECTOR method:
    x = 0:pi/4:(2*pi - pi/4);
    y = mean_r(n, :);
    [tuning_arr(n, :), res_arr(n, :)] = nlinfit(x, y, cos_fun, [1 1 0]);
    fit_arr(n, :) = cos_fun(tuning_arr(n, :), ang_fit); 
    [mx_arr, pref_idx] = max(fit_arr(n, :), [], 2);
    pref_arr(n) = ang_fit(pref_idx);
end
save ML mean_r sd_r
save PV tuning_arr pref_arr

% QUESTION 2 - plot histogram of preferred angles of all neurons by 
figure;
histogram(pref_arr, 'BinEdges', 0:pi/4:(2*pi));
title('Histogram of Preferred Angles for All 143 Neurons');
xlabel('Preferred Direction (in radians, 2pi = 6.28rad');
ylabel('Number of Neurons');

% QUESTION 2 - Fit of cosine curve
res2_tot_arr = zeros(N, 8);
for n = 1:N
    for dir = 1:8
        for trial = 1:trials
            if fr(trial, n, dir) > 0
                res2_tot_arr(n, dir) = res2_tot_arr(n, dir) + (fr(trial, n, dir) - mean_r(n, dir))^2;
            end
        end
    end
end
res2_arr = sum(res_arr.^2, 2);
% remove artifacts from non-firing neurons
res2_arr(res2_arr<1e-25) = NaN;
[res2_mn, res2_idx] = min(res2_arr);
r2_tot_arr = 1 - (res2_arr ./ sum(res2_tot_arr, 2));
[r2_mn, r2_idx] = max(r2_tot_arr);

% chi-squared test for preferred directions
e = N/8*ones(8, 1);
[h, p, stats] = chi2gof(pref_arr, 'Expected', e, 'Edges', 0.5:1:8.5);
fprintf('p-value of preferred directions: %d \n', p);

