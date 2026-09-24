function run_all()
%RUN_ALL  Every smoke test in this toolkit.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root,'core'), fullfile(root,'drivers'), here);
d = dir(fullfile(here, '*_smoke.m'));
fails = {}; t0 = tic;
for k = 1:numel(d)
    [~, n] = fileparts(d(k).name);
    t = tic;
    try
        evalc([n '()']);
        fprintf('PASS %-34s %5.1fs\n', n, toc(t));
    catch ME
        fprintf('FAIL %-34s %5.1fs  %s\n', n, toc(t), strtok(ME.message, newline));
        fails{end+1} = n; %#ok<AGROW>
    end
    close all force; drawnow;
end
fprintf('\n%d/%d passed in %.0f s\n', numel(d)-numel(fails), numel(d), toc(t0));
if ~isempty(fails), fprintf('failed: %s\n', strjoin(fails, ', ')); end
end
