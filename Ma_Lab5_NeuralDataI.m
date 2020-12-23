%% Neural Data Analysis I
% Student: Vanessa Ma
% import data
load('Lab5_CenterOutTrain.mat');
warning('off');

trials = length(go);
N = length(unit);

% create array for direction <> plot position
pos_arr = [6 3 2 1 4 7 8 9];
edges_psth = -1: .1 : 1;

% select neuron
rng(12456);
neuron = randi(N);

% plot rasters for chosen neuron
% figure
% hold on
% psth = zeros(8, length(edges_psth) - 1) ;
% for dir = 1:length(pos_arr)
%     subplot(3, 3, pos_arr(dir));
%     trial_idx = find(direction == dir);
%     for trial = 1:length(trial_idx)
%         cue = go(trial_idx(trial));
%         trial_arr = unit(neuron).times - cue;
%         trial_t = trial_arr(trial_arr < 1 & trial_arr > -1);
%         psth(dir, :) = psth(dir, :) + histcounts(trial_t, edges_psth);
%         for i = 1:length(trial_t)
%             line([trial_t(i) trial_t(i)], [1 2])
%             ylabel(sprintf('Direction %d', dir));
%         end        
%     end
% ylim([0 3])
% end
% subplot(3, 3, 2)
% title(sprintf('Raster Plot of 8 Directions for Neuron %d', neuron));
% subplot(3, 3, 8)
% xlabel('Time centered around go cue (s)');
% hold off

% plot rasters for chosen neuron 
figure 
hold on 
psth = zeros(8, length(edges_psth) - 1) ; 
for dir = 1:length(pos_arr)
    subplot(3, 3, pos_arr(dir)); 
    trial_idx = find(direction == dir); m=1; 
    for trial = 1:length(trial_idx) 
        cue = go(trial_idx(trial)); 
        trial_arr = unit(neuron).times - cue; 
        trial_t = trial_arr(trial_arr < 1 & trial_arr > -1); 
        psth(dir, :) = psth(dir, :) + histcounts(trial_t, edges_psth); 
        for i = 1:length(trial_t) 
            line([trial_t(i) trial_t(i)], [m-1 m]) 
            ylabel(sprintf('Direction %d', dir)); 
        end
        m=m+1;
    end
    ylim([0 length(trial_idx)]) 
end

figure;
for dir = 1:length(pos_arr)
    subplot(3, 3, pos_arr(dir));
    histogram('BinEdges', edges_psth, 'BinCounts', psth(dir, :));
    ylabel(sprintf('Direction %d', dir));
    ylim([0 40]);
    xlim([-1.1 1.1]);
end
subplot(3, 3, 2)
title(sprintf('Peri-Stimulus Time Plot of 8 Directions for Neuron: %d', neuron));
subplot(3, 3, 8)
xlabel('Time centered around go cue (s)');
hold off

% compute mean firing rate for each of your 143 neurons in each of the 8
% directions for a 2-second epoch centered on the go cue time.
% fit a cosine tuning curve for the mean firing rate as a function of the
% target angle.

% mean firing rates for all neurons in all directions
    % initialize 3D array of trials x directions x neurons
mean_r = zeros(N, 8);
sd_r = zeros(N, 8);
fr = zeros(trials, N, 8);
for n = 1:N
    if size(unit(n).times) > 0
        for dir = 1:8
            trial_idx = find(direction == dir);
            for trial = 1:length(trial_idx)
                cue = go(trial_idx(trial));
                trial_rel = unit(n).times - cue;
                trial_clean = trial_rel(trial_rel < 1 & trial_rel > -1);
                fr(trial_idx(trial), n, dir) = length(trial_clean) / 2;
            end
            mean_r(n, dir) = mean(nonzeros(fr(:, n, dir)));
            sd_r(n, dir) = std(nonzeros(fr(:, n, dir)));
        end
    end
end

% set up basic cosine function
cos_fun = @(p, theta) p(1) + p(2) *cos(theta - p(3));
    % set x as an array of directions (i.e. the target angle)
tuning_arr = zeros(N, 3);
fit_arr = zeros(N, 8);
res_arr = zeros(N, 8);
for n = 1:N
    x = 0:45:315;
    y = mean_r(n, :);
    [tuning_arr(n, :), res_arr(n, :)] = nlinfit(x, y, cos_fun, [1 1 0]);
    fit_arr(n, :) = cos_fun(tuning_arr(n, :), x); 
end

% determine best neuron
% res_arr(isnan(res_arr)) = 0;
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
fprintf('Best neuron by r2-value: %d \n', r2_idx);
fprintf('Corresponding r2 value: %d \n', r2_mn);
fprintf('Corresponding sum squared residuals value: %d \n \n', res2_arr(r2_idx));

fprintf('Best neuron by sum squared residuals: %d \n', res2_idx);
fprintf('Corresponding r2 value: %d \n', r2_tot_arr(res2_idx));
fprintf('Corresponding sum squared residuals value: %d \n', res2_mn);

% plot for best neuron by r2-value
figure;
plot(x, mean_r(r2_idx, :), 'b');
hold on
plot(x, fit_arr(r2_idx, :), 'k');
hold on

% find preferred angle for each neuron
[mx_arr, pref_arr] = max(fit_arr, [], 2);

% plot preferred direction for best neuron by r2-value
plot((pref_arr(r2_idx)-1)*45, mean_r(r2_idx, pref_arr(r2_idx)), 'kx', 'LineWidth', 2);
legend('Mean Firing Rates', 'Fit Points', 'Preferred Direction', 'location', 'northeast');
title(sprintf('Fit and Preferred Direction of Best Neuron by R2 %d', r2_idx));
xlabel('Degrees (0deg = direction 1)')
ylabel('Firing Rate (spikes / s)')
hold off

% plot preferred direction for best neuron by sum of squared errors
figure;
plot(x, mean_r(res2_idx, :), 'b');
hold on
plot(x, fit_arr(res2_idx, :), 'k');
hold on
plot((pref_arr(res2_idx)-1)*45, mean_r(res2_idx, pref_arr(res2_idx)), 'kx', 'LineWidth', 2);
legend('Mean Firing Rates', 'Fit Points', 'Preferred Direction', 'location', 'northwest ');
title(sprintf('Fit and Preferred Direction of Best Neuron by Sum Squared Errors %d', res2_idx));
xlabel('Degrees (0deg = direction 1)')
ylabel('Firing Rate (spikes / s)')
hold off

% plot histogram of preferred angles of all neurons
figure;
histogram(pref_arr, 'BinEdges', 0.5:1:8.5);
title('Histogram of Preferred Angles for All 143 Neurons');
xlabel('Directions');
ylabel('Number of Neurons');
