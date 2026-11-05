function data = build_and_simulate(cfg)
%BUILD_AND_SIMULATE Build a Simscape Multibody model of a manufacturer robot
% description and generate an identification dataset.
%
%   data = build_and_simulate(cfg)
%
% The model is created with smimport from the rigidBodyTree returned by
% panda_robot.m: the Franka Emika link geometry and meshes of loadrobot
% ('frankaEmikaPanda') with the masses, centre-of-mass positions and
% inertia tensors identified by Gaz et al. (2019). Every revolute joint is
% motion-actuated by the excitation trajectory (position, velocity and
% acceleration supplied analytically) and Simscape computes the actuator
% torques. Joint limits are disabled so that no limit-stop torques enter
% the data. The logged signals are converted to the MDH joint convention of
% the symbolic model:  q = S*q_m + delta,  tau = S*tau_m.
%
% smimport transfers only the mass properties (Inertia blocks). The link
% meshes are therefore exported to simulink/meshes and attached as massless
% File Solid blocks, so that the robot is visible in Mechanics Explorer
% without changing the dynamics. The model is organised as
%   Clock -> Excitation Trajectory -> <robot> Simscape Multibody -> Data Logging
% with the joint actuation and sensing grouped in their own subsystems.
%
% cfg fields:
%   robot        loadrobot name, e.g. 'universalUR5e'
%   model_name   name of the Simulink model to create
%   folder       folder in which the .slx is saved (meshes in folder/meshes)
%   label        robot name shown in the model, e.g. 'UR5e'
%   S, delta     joint-convention mapping (1xn)
%   q_min,q_max  joint limits in manufacturer coordinates (rad)
%   v_max        joint speed limits (rad/s)
%   frac         fraction of the half range and of v_max that may be used
%   f, phi       frequencies (Hz) and phases (rad) (1xn)
%   T_ramp       quintic ramp duration (s)
%   T_end, Ts    simulation length and sampling period (s)

    robot = panda_robot();      % loadrobot('frankaEmikaPanda') with published inertial data
    mdl = cfg.model_name;
    if ~isfield(cfg, 'label'), cfg.label = mdl; end
    if bdIsLoaded(mdl), close_system(mdl, 0); end
    meshdir = fullfile(cfg.folder, 'meshes');
    meshes = export_link_meshes(robot, meshdir, strrep(cfg.label, ' ', ''));
    addpath(meshdir);
    % build and simulate in a short temporary folder (Windows path limit for
    % Simulink build files); only the final .slx is saved to cfg.folder
    old = pwd;  cleanup = onCleanup(@() cd(old));
    work = fullfile(tempdir, 'sm_build');
    if ~isfolder(work), mkdir(work); end
    cd(work);
    target = fullfile(cfg.folder, [mdl '.slx']);
    if isfile(target), delete(target); end
    warning('off', 'all');
    smimport(robot, 'ModelName', mdl);
    warning('on', 'all');

    joints = find_system(mdl, 'SearchDepth', 1, 'ReferenceBlock', 'sm_lib/Joints/Revolute Joint');
    joints = sort_joints(joints, robot);
    nJ = numel(joints);
    add_visual_solids(mdl, meshes);

    % ---- trajectory block --------------------------------------------------
    [qc, dq] = excitation_amplitudes(cfg);
    traj = sprintf('%s/Excitation Trajectory', mdl);
    add_block('simulink/User-Defined Functions/MATLAB Function', traj, 'Position', [40 40 160 120]);
    code = trajectory_code(qc, dq, 2*pi*cfg.f, cfg.phi, cfg.T_ramp);
    rt = sfroot;  ch = rt.find('-isa', 'Stateflow.EMChart', 'Path', traj);
    ch.Script = code;
    add_block('simulink/Sources/Clock', [mdl '/Clock'], 'Position', [-40 70 -20 90]);
    add_line(mdl, 'Clock/1', 'Excitation Trajectory/1');

    % ---- reference bus: [q_d; qd_d; qdd_d] -> Goto [ref] ----------------------
    add_block('simulink/Signal Routing/Mux', [mdl '/Mux_ref'], 'Inputs', '3', 'Position', [200 40 205 120]);
    add_block('simulink/Signal Routing/Goto', [mdl '/Goto_ref'], 'GotoTag', 'ref', 'Position', [230 70 270 90]);
    for s = 1:3, add_line(mdl, sprintf('Excitation Trajectory/%d', s), sprintf('Mux_ref/%d', s)); end
    add_line(mdl, 'Mux_ref/1', 'Goto_ref/1');

    % ---- per-joint actuation and sensing: Joint i Interface -------------------------
    % From [ref] -> Joint i Interface (select q_i, qd_i, qdd_i; Simulink-PS converter
    % with user-provided derivatives; PS-Simulink converters for q, w, b, t)
    % -> Goto [y_i] = [q_i; qd_i; qdd_i; tau_i]
    units = {'rad', 'rad/s', 'rad/s^2', 'N*m'};
    sigs = {'q', 'qd', 'qdd', 'tau'};
    for j = 1:nJ
        blk = joints{j};  bn = get_param(blk, 'Name');
        set_param(blk, 'LowerLimitSpecify', 'off', 'UpperLimitSpecify', 'off', ...
            'MotionActuationMode', 'InputMotion', 'TorqueActuationMode', 'ComputedTorque', ...
            'SensePosition', 'on', 'SenseVelocity', 'on', 'SenseAcceleration', 'on', ...
            'SenseTorqueForce', 'on');
        ph = get_param(blk, 'PortHandles');
        x0 = 300 + 400*j;
        fr = sprintf('%s/From_ref_%d', mdl, j);
        add_block('simulink/Signal Routing/From', fr, 'GotoTag', 'ref', 'Position', [x0 400 x0+40 420]);
        grp = [];
        cin = sprintf('%s/In_%d', mdl, j);
        add_block('nesl_utility/Simulink-PS Converter', cin, 'Position', [x0+150 380 x0+180 440]);
        set_param(cin, 'Unit', 'rad', 'FilteringAndDerivatives', 'provide', 'UdotUserProvided', '2');
        cph = get_param(cin, 'PortHandles');
        for k = 1:3
            sel = sprintf('%s/Sel_%d_%s', mdl, j, sigs{k});
            add_block('simulink/Signal Routing/Selector', sel, 'NumberOfDimensions', '1', ...
                'IndexOptions', 'Index vector (dialog)', 'Indices', num2str((k-1)*nJ + j), ...
                'InputPortWidth', num2str(3*nJ), 'Position', [x0+80 370+25*k x0+110 385+25*k]);
            add_line(mdl, [get_param(fr, 'Name') '/1'], [get_param(sel, 'Name') '/1']);
            add_line(mdl, get_param(sel, 'PortHandles').Outport(1), cph.Inport(k));
            grp(end+1) = get_param(sel, 'Handle'); %#ok<AGROW>
        end
        add_line(mdl, cph.RConn(1), ph.LConn(end));
        grp(end+1) = get_param(cin, 'Handle'); %#ok<AGROW>
        mx = sprintf('%s/Mux_y_%d', mdl, j);
        add_block('simulink/Signal Routing/Mux', mx, 'Inputs', '4', 'Position', [x0+280 380 x0+285 460]);
        % sensing outputs: RConn(2:5) = q, w, b, t (verified in build_check)
        for s = 1:4
            cout = sprintf('%s/Out_%d_%s', mdl, j, sigs{s});
            add_block('nesl_utility/PS-Simulink Converter', cout, ...
                'Position', [x0+220 375+20*s x0+250 390+20*s]);
            set_param(cout, 'Unit', units{s});
            oph = get_param(cout, 'PortHandles');
            add_line(mdl, ph.RConn(1+s), oph.LConn(1));
            add_line(mdl, oph.Outport(1), get_param(mx, 'PortHandles').Inport(s));
            grp(end+1) = get_param(cout, 'Handle'); %#ok<AGROW>
        end
        grp(end+1) = get_param(mx, 'Handle'); %#ok<AGROW>
        gt = sprintf('%s/Goto_y_%d', mdl, j);
        add_block('simulink/Signal Routing/Goto', gt, 'GotoTag', sprintf('y%d', j), ...
            'Position', [x0+320 410 x0+360 430]);
        add_line(mdl, [get_param(mx, 'Name') '/1'], [get_param(gt, 'Name') '/1']);
        io = make_subsystem(mdl, grp, sprintf('Joint %d Interface', j));
        rename_ports(io, 'Inport', {'ref'});
        rename_ports(io, 'Outport', {sprintf('y%d', j)});
        name_connection_ports(io, j);
        fprintf('  joint %d -> %s\n', j, bn);
    end

    % ---- measurement bus: From [y_1..y_n] -> q, qd, qdd, tau -> To Workspace --
    add_block('simulink/Signal Routing/Mux', [mdl '/Mux_y'], 'Inputs', num2str(nJ), ...
        'Position', [3200 380 3205 560]);
    for j = 1:nJ
        fy = sprintf('%s/From_y_%d', mdl, j);
        add_block('simulink/Signal Routing/From', fy, 'GotoTag', sprintf('y%d', j), ...
            'Position', [3120 370+30*j 3160 390+30*j]);
        add_line(mdl, sprintf('From_y_%d/1', j), sprintf('Mux_y/%d', j));
    end
    for s = 1:4
        sel = sprintf('%s/Sel_%s', mdl, sigs{s});
        add_block('simulink/Signal Routing/Selector', sel, 'NumberOfDimensions', '1', ...
            'IndexOptions', 'Index vector (dialog)', 'Indices', sprintf('%d:4:%d', s, 4*nJ), ...
            'InputPortWidth', num2str(4*nJ), 'Position', [3260 340+60*s 3290 360+60*s]);
        add_line(mdl, 'Mux_y/1', sprintf('Sel_%s/1', sigs{s}));
        add_block('simulink/Sinks/To Workspace', sprintf('%s/Log_%s', mdl, sigs{s}), ...
            'VariableName', ['m_' sigs{s}], 'SaveFormat', 'Array', 'SampleTime', num2str(cfg.Ts), ...
            'Position', [3400 340+60*s 3480 360+60*s]);
        add_line(mdl, sprintf('Sel_%s/1', sigs{s}), sprintf('Log_%s/1', sigs{s}));
    end

    organize_model(mdl, cfg.label, robot, joints);
    set_param(mdl, 'StopTime', num2str(cfg.T_end), 'SolverType', 'Fixed-step', ...
        'Solver', 'ode3', 'FixedStep', num2str(cfg.Ts), ...
        'SimMechanicsOpenEditorOnUpdate', 'on', ...
        'InitFcn', 'addpath(fullfile(fileparts(get_param(bdroot, ''FileName'')), ''meshes''));');
    save_system(mdl, target);

    % ---- simulate and convert to the MDH convention -------------------------
    out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
    m.q = out.m_q; m.qd = out.m_qd; m.qdd = out.m_qdd; m.tau = out.m_tau;
    N = size(m.q, 1);
    data.time_data = (0:N-1).'*cfg.Ts;
    data.t_data   = m.q.*cfg.S + cfg.delta;
    data.td_data  = m.qd.*cfg.S;
    data.tdd_data = m.qdd.*cfg.S;
    data.T_data   = m.tau.*cfg.S;
    data.manufacturer = m;            % raw signals in manufacturer coordinates
    data.cfg = cfg;
    data.amplitude = dq;  data.center = qc;
    close_system(mdl, 0);
end

function add_visual_solids(mdl, meshes)
% Attach each link mesh as a massless File Solid to the body reference frame.
    palette = [0.24 0.55 0.37; 0.95 0.55 0.22; 0.18 0.53 0.76; 0.66 0.80 0.32;
               0.56 0.37 0.64; 0.55 0.58 0.60; 0.89 0.79 0.54];
    for k = 1:numel(meshes)
        sub = [mdl '/' meshes(k).body];
        if getSimulinkBlockHandle(sub) < 0, continue; end
        fs = [sub '/Visual'];
        if getSimulinkBlockHandle(fs) > 0
            % smimport already created the visual solid from the body's mesh
            % (detailed meshes given by file path); only the colour is set
            try
                set_param(fs, 'GraphicDiffuseColor', mat2str(palette(1 + mod(k-1, size(palette, 1)), :)));
            catch
            end
            continue;
        end
        add_block('sm_lib/Body Elements/File Solid', fs, 'Position', [420 200 480 260]);
        [~, name, ext] = fileparts(meshes(k).file);
        set_param(fs, 'ExtGeomFileName', [name ext], 'UnitType', 'Custom', 'ExtGeomFileUnits', 'm', ...
            'InertiaType', 'PointMass', 'Mass', '0', ...
            'GraphicType', 'FromGeometry', 'GraphicVisPropType', 'SimpleVisualProperties', ...
            'GraphicDiffuseColor', mat2str(palette(1 + mod(k-1, size(palette, 1)), :)), ...
            'GraphicOpacity', '1');
        ref = get_param([sub '/ReferenceFrame'], 'PortHandles');
        fph = get_param(fs, 'PortHandles');
        refport = [ref.LConn, ref.RConn];  solport = [fph.LConn, fph.RConn];
        add_line(sub, solport(1), refport(1), 'autorouting', 'on');
    end
end

function organize_model(mdl, label, robot, joints)
% Group the model into Excitation Trajectory -> robot -> Data Logging and lay
% out the robot subsystem: mechanism chain in one row, Joint i Interface blocks
% below their joints, reference and measurement buses at the sides.
    sigs = {'q', 'qd', 'qdd', 'tau'};
    jn = cellfun(@(b) get_param(b, 'Name'), joints, 'UniformOutput', false);
    keep = {'Clock', 'Excitation Trajectory', 'Log_q', 'Log_qd', 'Log_qdd', 'Log_tau'};
    blocks = find_system(mdl, 'SearchDepth', 1, 'Type', 'block');  blocks = blocks(2:end);
    rob = blocks(~ismember(get_param(blocks, 'Name'), keep));
    sub = make_subsystem(mdl, get_param(rob, 'Handle'), sprintf('%s (Simscape Multibody)', label));
    rename_ports(sub, 'Inport', {'q_d', 'qd_d', 'qdd_d'});
    rename_ports(sub, 'Outport', sigs);
    set_param(sub, 'BackgroundColor', 'orange');
    logs = make_subsystem(mdl, find_blocks(mdl, '^Log_'), 'Data Logging');
    rename_ports(logs, 'Inport', sigs);
    set_param(logs, 'BackgroundColor', 'lightBlue');
    set_param([mdl '/Excitation Trajectory'], 'BackgroundColor', 'cyan');
    try, Simulink.BlockDiagram.arrangeSystem(logs); catch, end

    % ---- top level ----------------------------------------------------------
    pos = {[-60 95 -30 125], [20 70 160 150], [240 50 460 170], [540 70 660 150]};
    blk = {'Clock', 'Excitation Trajectory', get_param(sub, 'Name'), 'Data Logging'};
    for i = 1:4, set_param([mdl '/' blk{i}], 'Position', pos{i}); end
    route_all(mdl);
    annotate(mdl, sprintf(['%s identification model. Every joint follows the excitation trajectory ' ...
        '(motion actuation), and Simscape Multibody computes the joint torques.\n' ...
        'Open the %s subsystem for the mechanism; the link meshes are massless visual solids, ' ...
        'the dynamics come from the published mass properties (Gaz et al., 2019).'], label, label), [-60 210 660 250]);

    % ---- robot subsystem ----------------------------------------------------
    chain = {robot.BaseName};
    for k = 1:robot.NumBodies
        chain = [chain, {matlab.lang.makeValidName(robot.Bodies{k}.Joint.Name), robot.Bodies{k}.Name}]; %#ok<AGROW>
    end
    chain = chain(cellfun(@(n) getSimulinkBlockHandle([sub '/' n]) > 0, chain));
    dx = 150;  y = 120;  xc = zeros(1, numel(chain));
    for i = 1:numel(chain)
        xc(i) = 260 + dx*(i-1);
        place([sub '/' chain{i}], xc(i), y);
    end
    % world, mechanism and solver configuration to the left of the base
    others = find_system(sub, 'SearchDepth', 1, 'Type', 'block');
    rb = get_param(others, 'ReferenceBlock');
    place(others{contains(rb, 'World Frame')}, 60, 60);
    place(others{contains(rb, 'Mechanism')}, 60, 140);
    place(others{contains(rb, 'Solver')}, 140, 210);
    % Joint i Interface below joint i
    for j = 1:numel(jn)
        x = xc(strcmp(chain, jn{j}));
        set_param(sprintf('%s/Joint %d Interface', sub, j), 'Position', [x-45 300 x+45 380], ...
            'BackgroundColor', 'lightBlue');
        set_param(sprintf('%s/From_ref_%d', sub, j), 'Position', [x-130 330 x-80 350], 'ShowName', 'off');
        set_param(sprintf('%s/Goto_y_%d', sub, j), 'Position', [x+80 330 x+130 350], 'ShowName', 'off');
    end
    % reference bus (left) and measurement bus (right)
    for s = 1:3
        set_param([sub '/' sigs{s} '_d'], 'Position', [20 260+50*s 50 274+50*s]);
    end
    set_param([sub '/Mux_ref'], 'Position', [120 300 125 420], 'ShowName', 'off');
    set_param([sub '/Goto_ref'], 'Position', [160 350 210 370], 'ShowName', 'off');
    xr = max(xc) + 220;
    for j = 1:numel(jn)
        set_param(sprintf('%s/From_y_%d', sub, j), 'Position', [xr 250+40*j xr+50 270+40*j], 'ShowName', 'off');
    end
    set_param([sub '/Mux_y'], 'Position', [xr+100 280 xr+105 530], 'ShowName', 'off');
    for s = 1:4
        set_param(sprintf('%s/Sel_%s', sub, sigs{s}), 'Position', [xr+170 250+60*s xr+210 270+60*s]);
        set_param([sub '/' sigs{s}], 'Position', [xr+280 253+60*s xr+310 267+60*s]);
    end
    route_all(sub);
    annotate(sub, sprintf('%s mechanism imported from the robot description (bodies, joints, mass properties, visual meshes)', label), [260 20 1200 40]);
    annotate(sub, 'Joint actuation and sensing: [ref] = [q_d; qd_d; qdd_d] drives each joint, [y_i] = [q_i; qd_i; qdd_i; tau_i] returns its sensed motion and torque', [260 420 1500 440]);
    for j = 1:numel(jn)
        io = sprintf('%s/Joint %d Interface', sub, j);
        try, Simulink.BlockDiagram.arrangeSystem(io); catch, end
    end
end

function name_connection_ports(io, j)
% name the physical ports of Joint j Interface after the joint ports they join
    map = {sprintf('In_%d', j), 'motion'; sprintf('Out_%d_q', j), 'q'; sprintf('Out_%d_qd', j), 'w';
           sprintf('Out_%d_qdd', j), 'b'; sprintf('Out_%d_tau', j), 't'};
    p = find_system(io, 'SearchDepth', 1, 'BlockType', 'PMIOPort');
    for i = 1:numel(p)
        pc = get_param(p{i}, 'PortConnectivity');
        nb = [pc.DstBlock, pc.SrcBlock];
        if isempty(nb), continue; end
        k = find(strcmp(map(:, 1), get_param(nb(1), 'Name')), 1);
        if ~isempty(k), set_param(p{i}, 'Name', map{k, 2}); end
    end
end

function place(blk, x, y)
    p = get_param(blk, 'Position');  w = p(3) - p(1);  h = p(4) - p(2);
    set_param(blk, 'Position', round([x - w/2, y - h/2, x + w/2, y + h/2]));
end

function route_all(sys)
    l = find_system(sys, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'line');
    try, Simulink.BlockDiagram.routeLine(l); catch, end
end

function annotate(sys, txt, pos)
    try
        a = Simulink.Annotation(sys, txt);  a.Position = pos;  a.FontSize = 11;
    catch
    end
end

function h = find_blocks(sys, expr)
    b = find_system(sys, 'SearchDepth', 1, 'Type', 'block');  b = b(2:end);
    h = get_param(b(~cellfun(@isempty, regexp(get_param(b, 'Name'), expr, 'once'))), 'Handle');
    h = [h{:}].';
end

function path = make_subsystem(parent, handles, name)
    if iscell(handles), handles = [handles{:}]; end
    before = find_system(parent, 'SearchDepth', 1, 'BlockType', 'SubSystem');
    Simulink.BlockDiagram.createSubsystem(handles(:).');
    after = find_system(parent, 'SearchDepth', 1, 'BlockType', 'SubSystem');
    new = setdiff(after, before);
    set_param(new{1}, 'Name', name);
    path = [parent '/' name];
end

function rename_ports(sub, type, names)
    p = find_system(sub, 'SearchDepth', 1, 'BlockType', type);
    n = cellfun(@(x) str2double(get_param(x, 'Port')), p);
    [~, o] = sort(n);  p = p(o);
    for i = 1:min(numel(p), numel(names))
        set_param(p{i}, 'Name', names{i});
    end
end

function [qc, A] = excitation_amplitudes(cfg)
% Amplitudes that keep |q - qc| <= frac*half range and |qd| <= frac*v_max
% over the whole trajectory, including the quintic ramp.
    qc = (cfg.q_max + cfg.q_min)/2;
    A  = cfg.frac*(cfg.q_max - cfg.q_min)/2;
    t  = (0:1e-3:cfg.T_end).';
    s  = min(t/cfg.T_ramp, 1);
    r  = 10*s.^3 - 15*s.^4 + 6*s.^5;
    rd = (30*s.^2 - 60*s.^3 + 30*s.^4)/cfg.T_ramp .* (t < cfg.T_ramp);
    w  = 2*pi*cfg.f;
    for j = 1:numel(A)
        unit_vel = max(abs(rd.*sin(w(j)*t + cfg.phi(j)) + r.*w(j).*cos(w(j)*t + cfg.phi(j))));
        A(j) = min(A(j), cfg.frac*cfg.v_max(j)/unit_vel);
    end
end

function joints = sort_joints(joints, robot)
% order the joint blocks along the kinematic chain
    names = {};
    for k = 1:robot.NumBodies
        jt = robot.Bodies{k}.Joint;
        if strcmp(jt.Type, 'revolute'), names{end+1} = jt.Name; end %#ok<AGROW>
    end
    bn = cellfun(@(b) get_param(b, 'Name'), joints, 'UniformOutput', false);
    idx = zeros(1, numel(names));
    for k = 1:numel(names)
        idx(k) = find(strcmp(bn, matlab.lang.makeValidName(names{k})) | strcmp(bn, names{k}), 1);
    end
    joints = joints(idx);
end

function code = trajectory_code(qc, dq, w, phi, Tr)
    v = @(x) mat2str(x(:).', 17);
    code = sprintf([ ...
        'function [q, qd, qdd] = fcn(t)\n' ...
        '%% Quintic-ramped multisine excitation with analytic derivatives\n' ...
        'qc = %s.'';  dq = %s.'';  w = %s.'';  phi = %s.'';  Tr = %s;\n' ...
        'if t <= 0, r = 0; rd = 0; rdd = 0;\n' ...
        'elseif t < Tr, s = t/Tr;\n' ...
        '    r = 10*s^3 - 15*s^4 + 6*s^5;\n' ...
        '    rd = (30*s^2 - 60*s^3 + 30*s^4)/Tr;\n' ...
        '    rdd = (60*s - 180*s^2 + 120*s^3)/Tr^2;\n' ...
        'else, r = 1; rd = 0; rdd = 0; end\n' ...
        'sn = sin(w*t + phi);  cs = cos(w*t + phi);\n' ...
        'q   = qc + dq.*r.*sn;\n' ...
        'qd  = dq.*(rd.*sn + r.*w.*cs);\n' ...
        'qdd = dq.*(rdd.*sn + 2*rd.*w.*cs - r.*w.^2.*sn);\n'], ...
        v(qc), v(dq), v(w), v(phi), mat2str(Tr, 17));
end
