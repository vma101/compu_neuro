% Stimulus Generation - Conjunction
function stim_plot_c(num, target, target_char, target_col, out_char, out_col)
% Output - create stimulus plot for conjunction search paradigm
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
    % (roughly even), divide up available 3 permutations other than target.
    % 2 methods of doing this are included to ensure randomization.
x_array=randperm(100, num);
y_array=randperm(100, num);
    % conjunction search, yes target condition
if(target == 1)
    if(num==4)
        stim_char = text(x_array(1), y_array(1), target_char);
        set(stim_char, 'color', target_col);
        stim_char = text(x_array(2), y_array(2), target_char);
        set(stim_char, 'color', out_col);
        stim_char = text(x_array(3), y_array(3), out_char);
        set(stim_char, 'color', out_col);
        stim_char = text(x_array(4), y_array(4), out_char);
        set(stim_char, 'color', target_col);
    else
        z = rand(1,1);
        if(z <= 0.5) % excess takes on target_char
            stim_char = text(x_array(1), y_array(1), target_char);
            set(stim_char, 'color', target_col);
            for i = 2:ceil(num/4)
                stim_char = text(x_array(i), y_array(i), target_char);
                set(stim_char, 'color', out_col);
            end
            for i = ceil(num/4 + 1):ceil(3 * num/4)
                stim_char = text(x_array(i), y_array(i), out_char);
                set(stim_char, 'color', target_col)
            end
            for i = ceil(3*num/4 + 1):num
                stim_char = text(x_array(i), y_array(i), out_char);
                set(stim_char, 'color', out_col);
            end
        else
            stim_char = text(x_array(1), y_array(1), out_char);
            set(stim_char, 'color', out_col);
            for i = 2:ceil(num/4)
                stim_char = text(x_array(i), y_array(i), out_char);
                set(stim_char, 'color', out_col);
            end
            for i = ceil(num/4 +1):ceil(3 * num/4)
                stim_char = text(x_array(i), y_array(i), target_char);
                set(stim_char, 'color', out_col)
            end
            for i = ceil(3*num/4 + 1):num
                stim_char = text(x_array(i), y_array(i), out_char);
                set(stim_char, 'color', target_col);
            end
        end
    end
    % conjunction search, no target condition
else
    if(num==4)
        for i = 1:2
            stim_char = text(x_array(i), y_array(i), target_char);
            set(stim_char, 'color', out_col);
        end
        stim_char = text(x_array(3), y_array(3), out_char);
        set(stim_char, 'color', target_col);
        stim_char = text(x_array(4), y_array(4), out_char);
        set(stim_char, 'color', out_col);
    else
        for i = 1:ceil(num/3)
            stim_char = text(x_array(i), y_array(i), out_char);
            set(stim_char, 'color', out_col);
        end
        for i = ceil(num/3 + 1):ceil(2 * num/3)
            stim_char = text(x_array(i), y_array(i), out_char);
            set(stim_char, 'color', target_col)
        end
        for i = ceil(2*num/3 + 1):num
            stim_char = text(x_array(i), y_array(i), target_char);
            set(stim_char, 'color', out_col);
        end
    end
end
end