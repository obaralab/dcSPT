function run_dcspt()
%RUN_DCSPT  Launch the dcSPT app — two colours, tracked independently, then compared.
%   Tabs: Cells | Detect | Co-motion | Pair.
%   Point it at a folder holding one cell's two colour stacks. It reads each file's own slice
%   labels to find out which pages are which colour, calibrates a detection threshold against each
%   colour's own noise, tracks them independently, and compares only CROSS-colour pairs.
%   The same shape as SPTinMatlab's launchers, so the two toolkits are driven the same way.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'core'), fullfile(here,'drivers'), fullfile(here,'app'));
dc_app();
end
