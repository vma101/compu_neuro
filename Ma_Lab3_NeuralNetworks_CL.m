%% Competitive Learning
% 2 categories: good and bad greebles
n_cat = 2;
% 3 features: boges, quiffs, dunths
n_features = 3;

data_good = xlsread('GoodGreeblesTraining.xls');
data_bad = xlsread('BadGreeblesTraining.xls');
data_training = [data_good; data_bad];
% target: 1/0 for 'Good' and 'Bad' Greebles
data_test = xlsread('Suspects_Test.xls');
figure;
plot_good_train = plot3(data_good(:, 1), data_good(:, 2), data_good(:, 3), 'ko');
hold on
plot_bad_train = plot3(data_bad(:, 1), data_bad(:, 2), data_bad(:, 3), 'go');
hold on

% random number generator is fixed - can comment out the following line to
% ensure randomness
rng(0, 'twister');
a = 0.01;
W = a.*rand(n_cat, n_features);
W_original = W;
quiver3(0, 0, 0, W(1, 1), W(1, 2), W(1, 3),  6, 'b-');
quiver3(0, 0, 0, W(2, 1), W(2, 2), W(2, 3),  6, 'r-');
disp('Original weights: ');
disp(W_original);

lr = 0.005;
run = 0;
W_diff = 6;
W_updates = 0.*W;
% re-run using training data as long as norm of differences in new and old
% weight matrices is not close to 0
while W_diff > 0.0000001
    W_old = W;
    run = run + 1;
    for i = 1:length(data_training)
        [W, idx] = train_cl(W, data_training(i, :)', lr);
        quiver3(0, 0, 0, W(1, 1), W(1, 2), W(1, 3),  6, 'r-');
        hold on
        quiver3(0, 0, 0, W(2, 1), W(2, 2), W(2, 3),  6, 'b-');
        hold on
    end
    W_diff = norm(W - W_old);
    % code to ensure that vector is not very very far off to address random
    % bias brought up during testing code in the lab instructions
    W_updates = W - W_old;
        % if the weight has not been updated for the whole 400-loop
        % iteration, then reset the weight
    for i = 1:size(W,1)
        for j = 1:size(W,2)
            if W(i, j) == W_old(i,j)
                W(i, j) = rand(1);
            end
        end
    end
end
    
disp('Final weights: ');
disp(W);

test_idx = test_cl(W, data_test');
data_test_result = [data_test'; test_idx];
% greeble classification for 'good' and 'bad' - 2/1 for bad/good
disp(data_test_result);
data_test_bad = data_test_result(: ,test_idx == 2)';
data_test_good = data_test_result(: , test_idx == 1)';

plot_good_test = plot3(data_test_good(:, 1), data_test_good(:, 2), data_test_good(:, 3), 'rx');
hold on
plot_bad_test = plot3(data_test_bad(:, 1), data_test_bad(:, 2), data_test_bad(:, 3), 'bx');
legend([plot_good_train plot_bad_train plot_good_test plot_bad_test], {'Good Greebles Training', 'Bad Greebles Training', 'Good Greebles Testing', 'Bad Greebles Testing'}, 'location', 'northeast');
hold off

figure;
quiver3(0, 0, 0, W_original(1, 1), W_original(1, 2), W_original(1, 3), 100, 'r-');
hold on
quiver3(0, 0, 0, W_original(2, 1), W_original(2, 2), W_original(2, 3), 100, 'b-');
hold on 
quiver3(0, 0, 0, W(1, 1), W(1, 2), W(1, 3), 'r--');
hold on
quiver3(0, 0, 0, W(2, 1), W(2, 2), W(2, 3), 'b--');
legend('Good Greeble Initial', 'Bad Greeble Initial', 'Good Greeble Final', 'Bad Greeble Final', 'location', 'northeast');
hold off

function [W_out, idx] = train_cl(W, input, lr)
% Competitive Learning Rule
% INPUTS
%   W - neural net, matrix object
%   input - array of input values for the feature neurons in W (column
%   vector taken transpose)
%   lr - learning rate parameter (between 0 and 1)
% OUTPUTS
%   W_out - updated weights for neural net
%   idx - array of classifications (1/2 for good/bad) for the input greeble
out = W*input;
[mx, idx] = max(out);
W(idx, :) = W(idx, :) + lr*(input' - W(idx, :));
W_out = W./repmat(sqrt(sum(W.^2, 2)), 1, size(W, 2));
end

function [idx] = test_cl(W, input)
% Competitive Learning Rule
% INPUTS
%   W - neural net, matrix object
%   input - array of input values for the feature neurons in W - whole
%   matrix is used, just transposed, because no adjustment is being made to
%   W
%   lr - learning rate parameter (between 0 and 1)
% OUTPUTS
%   idx - array of classifications (1/2 for good/bad) for the input
%   greebles
out = W*input;
[mx, idx] = max(out);
end
