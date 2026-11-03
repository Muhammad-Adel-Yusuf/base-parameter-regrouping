% Exact certificates for all twelve regroupings (rigid body / with actuation,
% symbolic / known dimensions, UR5e, IRB 120 and Panda). Writes exact_verification.txt
% into every results folder and a table to ../certificates_summary.txt.
here = fileparts(mfilename('fullpath'));  root = fileparts(here);
robots = {'UR5e', 'IRB120', 'Panda'};
cases = {'results', 'results/known_geometry', 'results/actuation', 'results/actuation_known_geometry'};
fid = fopen(fullfile(root, 'certificates_summary.txt'), 'w');
fprintf(fid, '%-7s %-34s %8s %10s %10s %9s  %s\n', 'robot', 'case', 'entries', 'nonzero', 'certified', 'time (s)', 'candidate');
for r = 1:numel(robots)
    addpath(fullfile(root, robots{r}));
    for c = 1:numel(cases)
        V = exact_verification(fullfile(root, robots{r}, cases{c}));
        fprintf(fid, '%-7s %-34s %8d %10d %10d %9.1f  %s\n', robots{r}, cases{c}, V.n_entries, ...
            V.n_nonzero_remainders, V.certified, V.run_time_s, V.candidate);
    end
    rmpath(fullfile(root, robots{r}));
end
fclose(fid);
type(fullfile(root, 'certificates_summary.txt'));
