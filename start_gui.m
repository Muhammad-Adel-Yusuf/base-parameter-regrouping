function app = start_gui()
%START_GUI Opens the interface for the certified base parameters of a serial manipulator.
%   From the project folder, type   start_gui   in the MATLAB Command Window.
    here = fileparts(mfilename('fullpath'));
    addpath(here, fullfile(here, 'app'), fullfile(here, 'mdh'), fullfile(here, 'core'), fullfile(here, 'robots'));
    app = base_parameter_app();
end
