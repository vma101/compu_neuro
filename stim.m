% Stimulus Generation
function [response, rt, validity, response_corr, cue_coord, target_coord] = stim(validity, box_array, ascii_array)
% INPUTS
% validity = 1/0/-1 for trial valid/invalid/neutral
% cue_position = tuple
% box_array = array of 4 rectangles
% ascii_array = array of ascii indexes for the right keypresses
% corresponding to each rectangle's position in box_array

% PROCESS
% 1) randomly generate cue position and target position (latter only for
% invalid trials)
% 2) set function variables such as cue time, isi time, stim time,
% response time 
% 3) get cue coordinates as output
% 4) use if statement to execute different displays depending on trial
% validity
% 5) all trials cue presentation code is the same, cue box takes on red
% outline for cue_time
% 6) for valid, invalid trials, target shows up after isi_time has lapsed,
% stays on for target_time, and target_coord records position of target.  
% After target disappears, respondents have response_time to press.
% getkeywait() function records reaction time (rt) and keypress (response).
% 7) check validity of response against ascii_array to return response_corr

% OUTPUT (array format)
% response - ascii keypress, left 0 if no button was pressed
% rt - reaction time (seconds)
% validity = 1/0/-1 for trial valid/invalid/neutral (same as input)
% response_corr - 1/0 for correct/incorrect
% cue_coord - tuple of x,y coordinates of the cue
% target_coord - tuple of x,y coordinates of the target

% choose box for cue and valid target
v_side = randi([1,4]);
i_side = ceil(rand*4);
response_corr = 0;
% choose box for invalid target
while(v_side == i_side)
    i_side = ceil(rand*4);
end

% other function variables
cue_time = 0.1;
isi_time = 0.2;
stim_time = 0.1;
response_time = 2;
target_coord = [0, 0];

% cue position return variable, same for all trials
cue_pos = get(box_array(v_side), 'Position');
cue_coord = [cue_pos(1) + 2, cue_pos(2) + 2];

%valid trial condition
if(validity == 1)
    box_array(v_side).LineWidth = 1;
    box_array(v_side).EdgeColor = 'r';
    pause(cue_time);
    box_array(v_side).LineWidth = 1;
    box_array(v_side).EdgeColor = 'k';
    pause(isi_time);
    target_pos = get(box_array(v_side), 'Position');
    target = text(target_pos(1) + 2, target_pos(2) + 2,'*', 'FontSize', 24);
    target_coord = [target_pos(1) + 2, target_pos(2) + 2];
    disp(target);
    pause(stim_time);
    delete(target);
    [response, rt] = getkeywait(response_time);
    if(response == ascii_array(v_side))
        response_corr = 1;
    end
% invalid trial condition
elseif(validity == 0)
    box_array(v_side).LineWidth = 1;
    box_array(v_side).EdgeColor = 'r';
    pause(cue_time);
    box_array(v_side).LineWidth = 1;
    box_array(v_side).EdgeColor = 'k';
    pause(isi_time);
    target_pos = get(box_array(i_side), 'Position');
    target = text(target_pos(1) + 2, target_pos(2) + 2,'*', 'FontSize', 24);
    target_coord = [target_pos(1) + 2, target_pos(2) + 2];
    disp(target);
    pause(stim_time);
    delete(target);
    [response, rt] = getkeywait(response_time);
    if(response == ascii_array(i_side))
        response_corr = 1;
    end
elseif(validity == -1)
    box_array(v_side).LineWidth = 1;
    box_array(v_side).EdgeColor = 'r';
    pause(cue_time);
    box_array(v_side).LineWidth = 1;
    box_array(v_side).EdgeColor = 'k';
    pause(isi_time);
    pause(stim_time);
    [response, rt] = getkeywait(response_time);
    if(response == -1)
        response_corr = 1;
    end
end
delete(box_array);
end
