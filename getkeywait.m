% Get a key within a time limit
function [CH, RT] = getkeywait(P)
%INPUTS
%   P - seconds to wait. If no key is pressed within P seconds, -1 is returned,
%   and if something went wrong during excution 0 is returned. Without 
%   argument, GETKEYWAIT waits until a key is pressed.
% OUTPUTS / PROCESS
%   CH - double representing the key pressed key as an ascii number, 
%   including backspace (8), space (32), enter (13), etc. If a non-ascii 
%   key (Ctrl, Alt, etc.) is pressed, CH will be NaN.  
%   RT - Reaction Time (s)

t00 = tic ; 
narginchk(0,1) ;

if nargin == 0
    P = -1 ;
elseif numel(P)~=1 || ~isnumeric(P) || ~isfinite(P) || P <= 0
    error('Argument should be a positive scalar.') ;
end

if P > 0
    % set up the timer
    tt = timer ;
    tt.timerfcn = 'uiresume' ;
    tt.startdelay = P ;
end

% Set up the figure
% May be the position property should be individually tweaked to avoid visibility
callstr = 'set(gcbf, ''Userdata'', double(get(gcbf,''Currentcharacter''))) ; uiresume ' ;
fh = figure(...
    'name', 'Press a key', ...
    'keypressfcn', callstr, ...
    'windowstyle', 'modal', ...
    'numbertitle', 'off', ...
    'position', [0 0  1 1], ... % really small in the corner
    'userdata', -1) ;

% start the timer, if needed
if P > 0
    % adjust the timer start delay
    tt.startdelay = max(0.001, round(P - toc(t00),3)) ;
    start(tt) ;
end

try
    % Wait for something to happen or the timer to run out    
    uiwait ;
    RT = toc(t00) ; % response time since start of this function
    CH = get(fh, 'userdata') ;
    if isempty(CH) % a non-ascii key was pressed, return a NaN
        CH = NaN ;
    end 
catch
    % something went wrong, return zero.
    CH = 0 ;
    RT = NaN ;
end

% clean up the timer ...
if P > 0
    stop(tt) ;
    delete(tt) ;
end

% ... and figure
if ishandle(fh)
    delete(fh) ;
end
end
