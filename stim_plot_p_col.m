% Stimulus Generation - Pop-Out (Color)
function stim_plot_p_col(num, target, target_char, target_col, out_char, out_col)
% Output - create stimulus plot for pop-out search paradigm
% Inputs
    % num = number of 'x', 'o' to plot
    % target = 1 / 0 for target presence
    % target_char = target character
    % target_col = target color
    % out_char = non-target character
    % out_col = non-target color
% Process
    % Create random arrays for x and y positions in the plots, range (0, 100)
    % Sort first by target presence, then by set size. 2 if-else statements
    % Set size difference - one condition for set size of 4, one condition
    % for anything larger than 4
    % In both conditions, if target present, set aside one marker as target
    % Set size = 4 - if target present, one of each color-char permutation
    % Set size > 4 - divide up set size into 3 portions 
    % (roughly even), divide up available 3 permutations other than target
x_array=randperm(100, num);
y_array=randperm(100, num);
    % pop-out search, yes target condition
if(target == 1)
    stim_char = text(x_array(1), y_array(1), target_char);
    set(stim_char, 'color', target_col);
    for i = 2:num
        stim_char = text(x_array(i), y_array(i), target_char);
        set(stim_char, 'color', out_col);
    end
    % pop-out search, no target condition
else
        for i = 1:num
            if(i <= num/2)
                char = target_char;
            else
                char = out_char;
            end
            stim_char = text(x_array(i), y_array(i), char);
            set(stim_char, 'color', out_col);
        end
end
end