function setup_path()
% SETUP_PATH  Puts engine and pols on the path, working from the location
% of this file rather than the current directory, so it does not matter
% which policy folder the session was started from.
%
% Call it once at the start of a session, before anything in engine or
% pols is used. The calibration files are deliberately left off the path:
% each policy folder carries its own copies and they are loaded from the
% working directory of whatever script is running.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'engine'));
addpath(fullfile(here, 'pols'));

end
