%% view_panda_simulation.m
% Opens the Panda Simscape Multibody model and plays the 40 s excitation
% trajectory in Mechanics Explorer at twice real time.
%
% Run it from the MATLAB desktop (not with -batch, -nodisplay or -nojvm):
% Mechanics Explorer needs Java and a display. If the window does not appear,
% press Ctrl+D in the model window (Update Diagram) or click the Mechanics
% Explorer button in the Simulink toolbar.

here = fileparts(mfilename('fullpath'));
cd(here);
mdl = 'Panda_identification';
open_system(mdl);
set_param(mdl, 'SimMechanicsOpenEditorOnUpdate', 'on', 'EnablePacing', 'on', 'PacingRate', 2);
sim(mdl);
