%% Supervised Learning Neural Networks - Back Propagation
% import kira data
k_1 = audioread('Kira_Training1.wav');
k_2 = audioread('Kira_Training2.wav');
k_3 = audioread('Kira_Training3.wav');
k_1_data = real(specgram(k_1));
k_2_data = real(specgram(k_2)); 
k_3_data = real(specgram(k_3));

% import pascal data
p_1 = audioread('Pascal_Training1.wav');
p_2 = audioread('Pascal_Training2.wav');
p_3 = audioread('Pascal_Training3.wav');
p_1_data = real(specgram(p_1));
p_2_data = real(specgram(p_2));
p_3_data = real(specgram(p_3));

% import test data
k_test = audioread('Test1.wav');
p_test = audioread('Test2.wav');
k_test_data = real(specgram(k_test));
p_test_data = real(specgram(p_test));

% define global variables for setup
n_group_time = size(k_1_data, 2);
n_chunks = round(size(k_1_data, 2)/n_group_time);
n_freqs = [46, 60; 61, 75; 76, 90; 91, 105; 106, 118; 119, 121];
n_inp = size(n_freqs, 1); %round(size(k_1_data, 1)/n_freq);
n_hidden = 4;
n_out = 1; % last layer should only have 1 node because it carries the decision
n_layers = 3;

% clean data
k_1_clean = clean_data(k_1_data, n_chunks, n_group_time, n_freqs, n_inp);
k_2_clean = clean_data(k_2_data, n_chunks, n_group_time, n_freqs, n_inp);
k_3_clean = clean_data(k_3_data, n_chunks, n_group_time, n_freqs, n_inp);
p_1_clean = clean_data(p_1_data, n_chunks, n_group_time, n_freqs, n_inp);
p_2_clean = clean_data(p_2_data, n_chunks, n_group_time, n_freqs, n_inp);
p_3_clean = clean_data(p_3_data, n_chunks, n_group_time, n_freqs, n_inp);
k_test_clean = clean_data(k_test_data, n_chunks, n_group_time, n_freqs, n_inp);
p_test_clean = clean_data(p_test_data, n_chunks, n_group_time, n_freqs, n_inp);

% collate all samples in terms of rows, to normalize along columns
% (features)
all_arr_raw = [k_1_clean'; k_2_clean'; k_3_clean'; p_1_clean'; p_2_clean'; p_3_clean'; k_test_clean'; p_test_clean'];
all_arr = normalize(all_arr_raw, 1);
% separate into sample / training and testing arrays
sample_arr = all_arr(1:6, :);
test_arr = all_arr(7:8, :);
% set targets for sample array
targets = [1 1 1 0 0 0];

% plots of sample features
    % first 3 features
figure;
plot3(sample_arr(1:3, 1), sample_arr(1:3, 2), sample_arr(1:3, 3), 'ro');
hold on
plot3(sample_arr(4:6, 1), sample_arr(4:6, 2), sample_arr(4:6, 3), 'bo');
title('First 3 features');
legend('Kira', 'Pascal', 'location', 'northwest');
hold off
    % latter 3 features
figure;
plot3(sample_arr(1:3, 4), sample_arr(1:3, 5), sample_arr(1:3, 6), 'r*');
hold on
plot3(sample_arr(4:6, 4), sample_arr(4:6, 5), sample_arr(4:6, 6), 'b*');
legend('Kira', 'Pascal', 'location', 'northwest');
title('Second 3 features');
hold off

%%% SETUP
% initialize weight matrices
rng(12345);
Wh = rand(n_hidden, n_inp); % goes from inp to hidden
Wo = rand(n_out, n_hidden); % goes from hidden to out
% initialize biases
bh = zeros(n_hidden, 1);
bo = zeros(n_out, 1);
% initialize parameters
lr = 0.05;
alpha = -1;
% collate Ws and bs into 3 separate arrays
W_arr = {0, Wh, Wo};
b_arr = {0, bh, bo};
% functions
    % sigmoid!
sigmoid_fn = @(x, alpha) 1./(1+exp(alpha*x));
dsigmoid_fn = @(x, alpha) sigmoid_fn(x, alpha).*(1 - sigmoid_fn(x, alpha));
    % analytic!
analytic_fn = @(x, alpha) log(1 + exp(x));
danalytic_fn = @(x, alpha) exp(x) ./ (1+ exp(x));
    % error function
error_fn = @(target, out) 0.5*(target - out).^2;
derror_fn = @(target, out) (target - out);
% choose transfer function
tFun = sigmoid_fn;
dFun = dsigmoid_fn;

% plots of samples prior to training
figure;
points = zeros(1, size(targets, 2));
for item = 1:size(targets, 2)
    inp = sample_arr(item, :)';
    [inp_arr, out_arr] = forward(inp, W_arr, b_arr, tFun, alpha);
    if item <= 3
        points(item) = plot(item, out_arr{n_layers}, 'r*');
    else
        points(item) = plot(item, out_arr{n_layers}, 'b*');
    end
    hold on
end
title('Plot of Samples Classification Before Training');
legend(points(3:4), {'Kira Sample', 'Pascal Sample'}, 'location', 'west');
xlabel('Sample - odd = Kira, even = Pascal');
ylabel('Classification - 1 = Kira, 0 = Pascal');
hold off    

%%% TRAINING
figure;
for iter = 1:10000
    % feed samples in one at a time through an array
    error_arr = zeros(1, length(targets));
    for item = 1:length(targets)
        target = targets(item);
        % feeding inputs as column arrays
        inp = sample_arr(item, :)';
        % FORWARD (starting layer 2)
        [inp_arr, out_arr] = forward(inp, W_arr, b_arr, tFun, alpha);
            % net input of layer appended to inp_arr, output of layer appended to
            % output_arr
            % reset inp to the output of the last layer for next layer
        out = out_arr{n_layers};
        error_tot = error_fn(target, out);
        % BACK PROPAGATION
        [W_arr, b_arr] = backward(W_arr, b_arr, inp_arr, out_arr, dFun, derror_fn, alpha, lr, target);
        error_arr(item) = error_tot;
        hold on
    end
    plot(iter, sum(error_arr), 'r*');
end
title('Error per iteration');
xlabel('Iteration');
ylabel('Halved sum squared error');
hold off
W_arr_final = W_arr;
b_arr_final = b_arr;

% plots of samples after training
figure;
points = zeros(1, size(targets, 2));
for item = 1:size(targets, 2)
    inp = sample_arr(item, :)';
    [inp_arr, out_arr] = forward(inp, W_arr, b_arr, tFun, alpha);
    if item <= 3
        points(item) = plot(item, out_arr{n_layers}, 'r*');
    else
        points(item) = plot(item, out_arr{n_layers}, 'b*');
    end
    hold on
end
legend(points(3:4), {'Kira Sample', 'Pascal Sample'}, 'location', 'west');
xlabel('Sample - odd = Kira, even = Pascal');
ylabel('Classification - 1 = Kira, 0 = Pascal');
title('Plot of Samples Classification After Training');
hold off   

%%% TESTING
figure;
k_inp = test_arr(1, :)';
[k_inp_arr, k_out_arr] = forward(k_inp, W_arr, b_arr, tFun, alpha);
k_out = k_out_arr{n_layers}; 
plot(1, k_out_arr{n_layers}, 'ro');
hold on
p_inp = test_arr(2, :)';
[p_inp_arr, p_out_arr] = forward(p_inp, W_arr, b_arr, tFun, alpha);
p_out = p_out_arr{n_layers}; 
plot(2, p_out_arr{n_layers}, 'bo');
legend('Kira Test', 'Pascal Test', 'location', 'northwest');
xlabel('Sample - 1 = Kira, 2 = Pascal');
ylabel('Classification - 1 = Kira, 0 = Pascal');
title('Plot of Test Samples Classification');
hold off
% calculate error on the tests (summed over the 2 samples)
test_error = error_fn(1, k_out) + error_fn(0, p_out);
disp(k_out);
disp(p_out);
disp(test_error);

function [inp_arr, out_arr] = forward(inp, W_arr, b_arr, tFun, alpha)
% INPUTS
    % inp = input into layer 1
    % W_arr = array of matrices for each layer
    % b_arr = array of biases for each layer
    % tFun = transfer function between net input and output of layer l
    % alpha = steepness of sigmoid
% PROCESS / OUTPUT
    % inp_arr = holds net input into each layer l (transformed by matrix and added bias)
    % out_arr = holds output of each layer l (apply transfer function)
inp_arr = cell(1, size(W_arr, 2));
out_arr = cell(1, size(W_arr, 2));
inp_arr{1} = inp;
out_arr{1} = inp;
for l = 2:size(W_arr, 2)
    inp_arr{l} = W_arr{l}*inp_arr{l-1} + b_arr{l};
    out_arr{l} = tFun(inp_arr{l}, alpha);
end
end

function [W_arr_new, b_arr_new] = backward(W_arr, b_arr, inp_arr, out_arr, dFun, derror_fn, alpha, lr, target)
% INPUTS (to implement for layer l)
    % W_arr = W matrix between layer 1 and layer l + 1
    % b_arr = array of biases
    % inp_arr = array of net inputs
    % out_arr = array of outputs 
    % dFun = transfer function's derivative
    % alpha = level of steepness for sigmoid
    % lr = learning rate
% PROCESS
    % 1) calculate dout/dinp (or here, dinp) for layer l, which is just the
    % transfer function's derivative
    % 2) calculate new delta for layer l (for use in layer l - 1), which is
    % 'error'
    % 3) update W by lr*delta*net_l-1 (i.e. the net input into layer l-1)
    % 4) update b by lr*delta    
% OUTPUT
    % W_arr_new - updated array of matrices
    % b_arr_new - updated array of biases
W_arr_new = W_arr;
b_arr_new = b_arr;
n_layers = size(W_arr, 2);
for layer = n_layers:-1:2
    dinp = dFun(inp_arr{layer}, alpha);
    if layer == n_layers
            % calculate derror and subsequent delta for update, next layer
        delta = derror_fn(target, out_arr{layer}) * dinp;
            % update Wo, bo and store in new array
    % for layers not top
    else
        % for given layer l, matrix to update is W_arr{l}, b to update is
        % b_arr{l}
        dinp = dFun(inp_arr{layer}, alpha);
        delta = (W_arr{layer + 1}' * delta).* dinp;
    end
        % update W, b for each layer and store in new array
    W_arr_new{layer} = W_arr{layer} + lr * delta * (out_arr{layer - 1})';
    b_arr_new{layer} = b_arr{layer} + lr * delta;
end
end

function data_clean = clean_data(data, cols, col_size, row_arr, inps)
% INPUTS
    % data = data to be cleaned
    % cols = how many cols to wind up with in the end
    % col_size = how many columns to average over in the original dataset
    % to wind up in the output data
    % row_arr = beginning and end rows to average across, in an array of
    % rows
    % inp_size = how many rows to wind up with in the end
data_clean = zeros(inps, cols);
[row_end, col_end] = size(data);
% clean columns (just averaging across)
data_col_clean = zeros(row_end, cols);
for i = 1:cols
    for j = 1:row_end
        data_col_clean(j, i) = mean(data(j, ((i - 1)*col_size + 1):min(i*col_size, col_end)));
    end
end
% clean rows (just averaging across)
for j = 1:inps
    for i = 1:cols
        data_clean(j, i) = mean(data_col_clean(row_arr(j, 1):min(row_arr(j, 2), row_end), i));
    end
end
end