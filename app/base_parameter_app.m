function app = base_parameter_app(varargin)
%BASE_PARAMETER_APP Interface for the certified base parameters of a serial manipulator.
%
%   base_parameter_app                  % opens the window
%   app = base_parameter_app('Visible','off')   % returns a handle struct (used by the tests)
%
% Tab "Robot":      number of joints, joint types (revolute / prismatic), the modified-DH table
%                   [alpha(i-1) a(i-1) d(i) theta0(i)], link lengths (known / unknown), actuation terms
%                   (rotor inertia, viscous and Coulomb friction) and the options of the run.
% Tab "Results":    summary of the certificates and the base parameters of each variant.
% Tab "Values":     standard (CAD) parameters of each link, the resulting numerical base parameters, and an
%                   independent Newton-Euler check of the certified regressor (Lagrange vs. recursive Newton-Euler).
% Tab "Identify":   least-squares identification of the base parameters, from loaded measured data or from a
%                   trajectory generated and simulated in the app (independent Newton-Euler reference, with
%                   optional Gaussian noise and outliers).
%
% Table entries are numbers or expressions: 0, pi/2, -pi/2, atan(3/4), 0.425, L1, L2+L3. Every name that is
% not a number is taken as a link length; "Detect lengths" lists them. Unknown lengths stay symbolic.
% A prismatic joint adds its displacement q_i to the offset d(i); a revolute joint adds q_i to theta0(i).

    ip = inputParser;  addParameter(ip, 'Visible', 'on');  parse(ip, varargin{:});
    here = fileparts(mfilename('fullpath'));  root = fileparts(here);
    addpath(fullfile(root, 'mdh'), fullfile(root, 'core'), fullfile(root, 'robots'));
    warning('off', 'mdh_table:alpha');
    folder_of = struct('symbolic', 'symbolic', 'known', 'known_geometry', ...
        'actuation', 'actuation', 'actuation_known', 'actuation_known_geometry');
    stdcols = {'m (kg)', 'Lcx', 'Lcy', 'Lcz', 'Ixx', 'Iyy', 'Izz', 'Ixy', 'Ixz', 'Iyz', 'Ia', 'Fv', 'Fc'};
    S = struct('mdh', [], 'out', [], 'variants', {{}}, 'message', '', 'ident', [], 'nominal', [], 'theta', struct());

    %% ---------------------------------------------------------------- window and layout
    fig = uifigure('Name', 'SYCIRO: from the MDH table to certified base parameters and identification', 'Position', [60 40 1240 780], ...
        'Visible', ip.Results.Visible);
    gfig = uigridlayout(fig, [2 1]);  gfig.Padding = [6 6 6 6];  gfig.RowHeight = {30, '1x'};
    banner = uilabel(gfig, 'FontWeight', 'bold', 'FontSize', 16, 'Text', ...
        'SYCIRO  |  SYmbolic Certified Identification of RObots:  MDH table  >  certified base parameters  >  identification');
    banner.Layout.Row = 1;
    tabs = uitabgroup(gfig);  tabs.Layout.Row = 2;
    tabRobot = uitab(tabs, 'Title', 'Robot');
    tabRes   = uitab(tabs, 'Title', 'Results');
    tabVal   = uitab(tabs, 'Title', 'Values');
    tabId    = uitab(tabs, 'Title', 'Identify');

    %% ---- Robot tab
    gR = uigridlayout(tabRobot, [1 2]);  gR.ColumnWidth = {380, '1x'};
    pl = uipanel(gR, 'Title', 'Robot and options');
    gl = uigridlayout(pl, [18 2]);  gl.RowHeight = repmat({26}, 1, 18);  gl.ColumnWidth = {150, '1x'};  gl.RowSpacing = 6;
    r = 0;
    r = r + 1;  lbl(gl, 'Name', r);
    edName = uieditfield(gl, 'text', 'Value', 'MyRobot');  edName.Layout.Row = r;  edName.Layout.Column = 2;
    r = r + 1;  lbl(gl, 'Number of joints', r);
    spN = uispinner(gl, 'Limits', [1 12], 'Value', 3, 'RoundFractionalValues', 'on', 'ValueChangedFcn', @(~, ~) onNChanged());
    spN.Layout.Row = r;  spN.Layout.Column = 2;
    r = r + 1;  lbl(gl, 'Load an example', r);
    ddPreset = uidropdown(gl, 'Items', presetList());  ddPreset.Layout.Row = r;  ddPreset.Layout.Column = 2;
    r = r + 1;
    btnLoad = uibutton(gl, 'Text', 'Load example into the tables', 'ButtonPushedFcn', @(~, ~) loadPreset(ddPreset.Value));
    btnLoad.Layout.Row = r;  btnLoad.Layout.Column = [1 2];
    r = r + 1;  lbl(gl, 'Link lengths', r);
    ddLen = uidropdown(gl, 'Items', {'Known (numeric values)', 'Unknown (symbolic)', 'Both'}, 'Value', 'Both');
    ddLen.Layout.Row = r;  ddLen.Layout.Column = 2;
    r = r + 1;  lbl(gl, 'Rotor inertia + friction', r);
    ddAct = uidropdown(gl, 'Items', {'Rigid body only', 'With rotor inertia and friction', 'Both'}, 'Value', 'Rigid body only');
    ddAct.Layout.Row = r;  ddAct.Layout.Column = 2;
    r = r + 1;  lblCompute = lbl(gl, 'What to compute', r);
    lblCompute.Tooltip = 'Base parameters: the full pipeline (regrouping and certificates). Dynamic model only: the symbolic torque equations, not regrouped.';
    ddCompute = uidropdown(gl, 'Items', {'Base parameters and certificates', 'Dynamic model only (symbolic equations)'}, ...
        'ValueChangedFcn', @(~, ~) onComputeChanged());
    ddCompute.Layout.Row = r;  ddCompute.Layout.Column = 2;
    r = r + 1;  lbl(gl, 'Derivation', r);
    ddDeriv = uidropdown(gl, 'Items', {'Euler-Lagrange', 'Newton-Euler'});
    ddDeriv.Layout.Row = r;  ddDeriv.Layout.Column = 2;
    ddDeriv.Tooltip = 'Dynamic model only: the method used to derive the equations. The base-parameter pipeline always uses Euler-Lagrange.';
    r = r + 1;
    cbPar = uicheckbox(gl, 'Text', 'Use a parallel pool (faster)', 'Value', license('test', 'Distrib_Computing_Toolbox') == 1);
    cbPar.Layout.Row = r;  cbPar.Layout.Column = [1 2];
    r = r + 1;
    cbRec = uicheckbox(gl, 'Text', 'Exact reconstruction check (slow above 6 joints)', 'Value', true);
    cbRec.Layout.Row = r;  cbRec.Layout.Column = [1 2];
    r = r + 1;
    cbMin = uicheckbox(gl, 'Text', 'Exact minimality proof', 'Value', true);
    cbMin.Layout.Row = r;  cbMin.Layout.Column = [1 2];
    r = r + 1;  lbl(gl, 'Results folder', r);
    edFolder = uieditfield(gl, 'text', 'Value', fullfile(root, 'results_app'));  edFolder.Layout.Row = r;  edFolder.Layout.Column = 2;
    r = r + 1;
    btnBrowse = uibutton(gl, 'Text', 'Choose results folder ...', 'ButtonPushedFcn', @(~, ~) onBrowse());
    btnBrowse.Layout.Row = r;  btnBrowse.Layout.Column = [1 2];
    r = r + 1;
    btnRun = uibutton(gl, 'Text', 'RUN: find and certify the base parameters', 'FontWeight', 'bold', ...
        'BackgroundColor', [0.80 0.92 0.80], 'ButtonPushedFcn', @(~, ~) runPipeline());
    btnRun.Layout.Row = r;  btnRun.Layout.Column = [1 2];
    r = r + 1;
    btnSave = uibutton(gl, 'Text', 'Save this robot as robots/mdh_<name>.m', 'ButtonPushedFcn', @(~, ~) saveRobotFile());
    btnSave.Layout.Row = r;  btnSave.Layout.Column = [1 2];
    r = r + 1;
    lblStatus = uilabel(gl, 'Text', 'Ready.', 'WordWrap', 'on');
    lblStatus.Layout.Row = [r r + 1];  lblStatus.Layout.Column = [1 2];

    pr = uigridlayout(gR, [4 1]);  pr.RowHeight = {22, '1x', 22, 190};
    uilabel(pr, 'Text', 'Modified DH table (Khalil-Kleinfinger). Entries: numbers or expressions such as pi/2, atan(3/4), L1, L2+L3.', ...
        'WordWrap', 'on');
    tblMDH = uitable(pr, 'ColumnName', {'Joint type', 'alpha(i-1)  [rad]', 'a(i-1)  [m]', 'd(i)  [m]', 'theta0(i)  [rad]'}, ...
        'ColumnEditable', true(1, 5), 'ColumnFormat', {{'Revolute', 'Prismatic'}, 'char', 'char', 'char', 'char'}, ...
        'RowName', 'numbered', 'ColumnWidth', {110, 140, 140, 140, 140});
    gd = uigridlayout(pr, [1 4]);  gd.ColumnWidth = {'1x', 130, 130, 130};  gd.Padding = [0 0 0 0];
    lblDim = uilabel(gd, 'Text', 'Link lengths (m); NaN = unknown', 'WordWrap', 'on');
    lblDim.Tooltip = 'Numeric values are needed only for the variants with known lengths.';
    uibutton(gd, 'Text', 'Detect lengths', 'ButtonPushedFcn', @(~, ~) detectDimensions());
    uibutton(gd, 'Text', 'Add row', 'ButtonPushedFcn', @(~, ~) addDimRow());
    uibutton(gd, 'Text', 'Remove last row', 'ButtonPushedFcn', @(~, ~) removeDimRow());
    tblDim = uitable(pr, 'ColumnName', {'Name', 'Value (m)'}, 'ColumnEditable', [true true], ...
        'ColumnFormat', {'char', 'numeric'}, 'RowName', 'numbered');

    %% ---- Results tab
    gS = uigridlayout(tabRes, [6 1]);  gS.RowHeight = {170, 40, 30, '1x', 24, 140};
    tblSum = uitable(gS, 'ColumnName', {'Variant', 'p', 'p_b', 'Dependency', 'Reconstruction', 'Minimality', ...
        'e_max (N m)', 'r_geom', 'r_RNEA', 'Time (s)'}, 'RowName', {});
    uilabel(gS, 'WordWrap', 'on', 'FontAngle', 'italic', 'Text', ...
        ['The certificates cover the algebraic reduction of the model defined by the entered table (dependency identity, ' ...
        'reconstruction of the symbolic torques, minimality). They do not verify the table, the rigid-body assumptions, or the real robot.']);
    gb = uigridlayout(gS, [1 5]);  gb.ColumnWidth = {80, 220, 150, 150, '1x'};  gb.Padding = [0 0 0 0];
    uilabel(gb, 'Text', 'Variant:');
    ddVar = uidropdown(gb, 'Items', {'-'}, 'ValueChangedFcn', @(~, ~) showVariant());
    uibutton(gb, 'Text', 'Open results folder', 'ButtonPushedFcn', @(~, ~) openFolder());
    uibutton(gb, 'Text', 'Save list as text', 'ButtonPushedFcn', @(~, ~) saveList());
    lblCount = uilabel(gb, 'Text', '');
    taParams = uitextarea(gS, 'Editable', 'off', 'FontName', 'Consolas');
    uilabel(gS, 'Text', 'Run log');
    taLog = uitextarea(gS, 'Editable', 'off', 'FontName', 'Consolas');

    %% ---- Values tab
    gV = uigridlayout(tabVal, [6 1]);  gV.RowHeight = {44, 190, 30, 24, '1x', 46};
    uilabel(gV, 'Text', ['Standard parameters of each link (optional): mass, centre of mass in the link frame, centroidal inertia ' ...
        'in the link frame, and for each joint the rotor inertia Ia and the viscous (Fv) and Coulomb (Fc) friction. ' ...
        'Select a variant in the Results tab, then compute the numerical base parameters.'], 'WordWrap', 'on');
    tblStd = uitable(gV, 'ColumnName', stdcols, 'ColumnEditable', true, 'RowName', 'numbered');
    gv2 = uigridlayout(gV, [1 5]);  gv2.ColumnWidth = {260, 180, 240, '1x', 10};  gv2.Padding = [0 0 0 0];
    uibutton(gv2, 'Text', 'Compute numerical base parameters', 'ButtonPushedFcn', @(~, ~) computeNominal());
    uibutton(gv2, 'Text', 'Save as CSV', 'ButtonPushedFcn', @(~, ~) saveTable(tblNom, 'base_parameter_values.csv'));
    uibutton(gv2, 'Text', 'Check against Newton-Euler', 'ButtonPushedFcn', @(~, ~) rneaCheck());
    lblNom = uilabel(gV, 'Text', '');
    tblNom = uitable(gV, 'ColumnName', {'#', 'Base parameter', 'Value'}, 'RowName', {}, 'ColumnWidth', {50, 700, 160});
    lblRnea = uilabel(gV, 'Text', ['Independent check: evaluates the certified Lagrange regressor and the recursive Newton-Euler ' ...
        'algorithm at the same random states and parameter values, and reports their largest difference.'], 'WordWrap', 'on');

    %% ---- Identify tab
    gI = uigridlayout(tabId, [6 1]);  gI.RowHeight = {60, 54, 30, 24, '1x', 24};
    uilabel(gI, 'Text', ['Least-squares identification of the base parameters of the selected variant (lengths known). ' ...
        'Data: a .mat file with N-by-n arrays q, qd, qdd, tau, or a .csv file with 4n columns [q qd qdd tau]; or generate an ' ...
        'excitation trajectory below, with torques from the independent Newton-Euler reference and optional noise and outliers. ' ...
        'Use the MDH convention of the table (joint force for a prismatic joint).'], 'WordWrap', 'on');
    gTraj = uigridlayout(gI, [1 6]);  gTraj.ColumnWidth = {150, 150, 150, 150, 150, '1x'};  gTraj.Padding = [0 0 0 0];
    [~, edDur]   = fieldBox(gTraj, 'Duration (s)', 20, [1 Inf]);
    [~, edFs]    = fieldBox(gTraj, 'Sample rate (Hz)', 50, [1 Inf]);
    [~, edAmp]   = fieldBox(gTraj, 'Amplitude (rad or m)', 0.6, [0 Inf]);
    [~, edNoise] = fieldBox(gTraj, 'Noise (% of torque std)', 2, [0 100]);
    [~, edOut]   = fieldBox(gTraj, 'Outliers (%)', 0, [0 100]);
    uibutton(gTraj, 'Text', 'Generate trajectory and identify', 'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(~, ~) generateAndIdentify());
    gi2 = uigridlayout(gI, [1 4]);  gi2.ColumnWidth = {200, 200, 180, '1x'};  gi2.Padding = [0 0 0 0];
    uibutton(gi2, 'Text', 'Load data and identify ...', 'ButtonPushedFcn', @(~, ~) onLoadData());
    uibutton(gi2, 'Text', 'Save estimates as CSV', 'ButtonPushedFcn', @(~, ~) saveTable(tblId, 'base_parameter_estimates.csv'));
    lblId = uilabel(gI, 'Text', '');
    tblId = uitable(gI, 'ColumnName', {'#', 'Base parameter', 'Estimate'}, 'RowName', {}, 'ColumnWidth', {50, 700, 160});
    lblId2 = uilabel(gI, 'Text', '');

    setN(3);
    onComputeChanged();
    setStatus('Ready. Enter the table, or load an example, then press RUN.');

    %% ---------------------------------------------------------------- returned handle struct
    app = struct('fig', fig, 'setN', @setN, 'loadPreset', @loadPreset, 'run', @runPipeline, ...
        'setOption', @setOption, 'setTable', @setTable, 'setDims', @setDims, 'setStandard', @setStandard, ...
        'computeNominal', @computeNominal, 'identify', @identify, 'detect', @detectDimensions, ...
        'saveRobot', @saveRobotFile, 'state', @getState, 'message', @getMessage, ...
        'rneaCheck', @rneaCheck, 'setTrajectory', @setTrajectory, 'generateAndIdentify', @generateAndIdentify);

    %% ================================================================ callbacks and helpers (nested)
    function items = presetList()
        d = dir(fullfile(root, 'robots', 'mdh_*.m'));
        nm = cellfun(@(s) s(5:end-2), {d.name}, 'UniformOutput', false);
        items = [{'(choose)'}, nm];
    end

    function s = dispExpr(e)
        % display form: decimals with 6 significant digits (the saved files keep the exact fractions)
        s = char(vpa(e, 6));
        s = regexprep(s, '(?<![\d.])1\.0\*', '');          % no redundant unit coefficient
    end

    function s = getState()
        s = S;                          % a nested function sees the live workspace (an anonymous function would copy S)
    end

    function m = getMessage()
        m = S.message;
    end

    function L = lbl(parent, txt, row)
        L = uilabel(parent, 'Text', txt);
        L.Layout.Row = row;  L.Layout.Column = 1;
    end

    function [L, F] = fieldBox(parent, txt, val, lims)
        % Compact label-over-field control for a row of numeric settings.
        p = uigridlayout(parent, [2 1]);  p.Padding = [0 0 0 0];  p.RowHeight = {16, '1x'};  p.RowSpacing = 2;
        L = uilabel(p, 'Text', txt, 'FontSize', 10);
        F = uieditfield(p, 'numeric', 'Value', val, 'Limits', lims);
    end

    function setStatus(msg)
        S.message = msg;
        if isvalid(lblStatus), lblStatus.Text = msg; end
    end

    function notify(msg)
        setStatus(msg);
        if strcmp(fig.Visible, 'on'), uialert(fig, msg, 'Please check', 'Icon', 'warning'); end
    end

    function setN(n)
        spN.Value = n;  onNChanged();
    end

    function onNChanged()
        n = round(spN.Value);
        D = tblMDH.Data;
        if isempty(D), D = cell(0, 5); end
        if size(D, 1) > n
            D = D(1:n, :);
        else
            for i = size(D, 1) + 1:n, D(i, :) = {'Revolute', '0', '0', '0', '0'}; end
        end
        tblMDH.Data = D;
        cbRec.Value = n <= 6;
        S_std = tblStd.Data;
        if isempty(S_std), S_std = zeros(0, numel(stdcols)); end
        if size(S_std, 1) > n, S_std = S_std(1:n, :); else, S_std(end+1:n, :) = 0; end
        tblStd.Data = S_std;
        if n >= 7, setStatus(sprintf('%d joints: symbolic extraction grows quickly; expect many minutes (7 joints: about 10 minutes with a parallel pool).', n)); end
    end

    function onBrowse()
        d = uigetdir(edFolder.Value, 'Results folder');
        if ischar(d) && d ~= 0, edFolder.Value = d; end
    end

    function setOption(name, value)
        switch lower(name)
            case 'name',           edName.Value = value;
            case 'lengths',        ddLen.Value = value;
            case 'actuation',      ddAct.Value = value;
            case 'parallel',       cbPar.Value = logical(value);
            case 'reconstruction', cbRec.Value = logical(value);
            case 'minimality',     cbMin.Value = logical(value);
            case 'folder',         edFolder.Value = value;
            case 'compute',        ddCompute.Value = value;  onComputeChanged();
            case 'derivation',     ddDeriv.Value = value;
            otherwise, error('Unknown option %s.', name);
        end
    end

    function setTable(D)
        % D: n-by-5 cell {type, alpha, a, d, theta0}; type 'Revolute' or 'Prismatic'
        spN.Value = size(D, 1);  tblMDH.Data = D;  onNChanged();
    end

    function setDims(names, values)
        tblDim.Data = [names(:), num2cell(values(:))];
    end

    function setStandard(M)
        tblStd.Data = M;
    end

    function loadPreset(name)
        if strcmp(name, '(choose)'), return; end
        try
            mdh = feval(['mdh_' name]);
        catch err
            notify(sprintf('Could not load %s: %s', name, err.message));  return
        end
        n = mdh.n;  D = cell(n, 5);
        for i = 1:n
            D{i, 1} = ternary(mdh.prismatic(i), 'Prismatic', 'Revolute');
            D{i, 2} = char(mdh.alpha(i));  D{i, 3} = char(mdh.a(i));
            D{i, 4} = char(mdh.d(i));      D{i, 5} = char(mdh.theta0(i));
        end
        edName.Value = mdh.name;
        setTable(D);
        tblDim.Data = [mdh.dim_names(:), num2cell(mdh.dim_values(:))];
        tblStd.Data = zeros(n, numel(stdcols));
        nomfun = ['nominal_parameters_' lower(name)];
        if exist(nomfun, 'file') == 2
            try
                [~, ~, c] = feval(nomfun);
                tblStd.Data = [c.m(:), c.Lc, c.Ic(:, 1:6), zeros(n, 3)];
            catch
            end
        end
        setStatus(sprintf('Example "%s" loaded: %d joints.', name, n));
    end

    function [rows, vars] = parseTable()
        D = tblMDH.Data;  n = size(D, 1);
        rows = sym(zeros(n, 5));
        for i = 1:n
            for c = 2:5
                txt = strtrim(char(string(D{i, c})));
                if isempty(txt), txt = '0'; end
                try
                    rows(i, c - 1) = str2sym(txt);
                catch
                    error('Row %d, column "%s": cannot read "%s". Use numbers or expressions such as pi/2, 0.425 or L1.', ...
                        i, tblMDH.ColumnName{c}, txt);
                end
            end
            rows(i, 5) = double(strcmp(D{i, 1}, 'Prismatic'));
        end
        vars = symvar(rows(:, 1:4));
    end

    function detectDimensions()
        try
            [~, vars] = parseTable();
        catch err
            notify(err.message);  return
        end
        names = arrayfun(@char, vars, 'UniformOutput', false);
        old = tblDim.Data;  vals = nan(numel(names), 1);
        for k = 1:numel(names)
            if ~isempty(old)
                j = find(strcmp(old(:, 1), names{k}), 1);
                if ~isempty(j) && ~isempty(old{j, 2}), vals(k) = old{j, 2}; end
            end
        end
        tblDim.Data = [names(:), num2cell(vals)];
        setStatus(sprintf('%d link length(s) found in the table: %s', numel(names), strjoin(names, ', ')));
    end

    function addDimRow()
        D = tblDim.Data;  if isempty(D), D = cell(0, 2); end
        D(end + 1, :) = {'L', NaN};  tblDim.Data = D;
    end

    function removeDimRow()
        D = tblDim.Data;  if ~isempty(D), tblDim.Data = D(1:end - 1, :); end
    end

    function mdh = buildMdh(needKnown)
        [rows, vars] = parseTable();
        names = arrayfun(@char, vars, 'UniformOutput', false);
        D = tblDim.Data;  vals = nan(1, numel(names));
        for k = 1:numel(names)
            if ~isempty(D)
                j = find(strcmp(D(:, 1), names{k}), 1);
                if ~isempty(j) && ~isempty(D{j, 2}) && isfinite(D{j, 2}), vals(k) = D{j, 2}; end
            end
        end
        if needKnown && any(isnan(vals))
            error('Enter numeric values for the link lengths: %s (or choose "Unknown (symbolic)").', strjoin(names(isnan(vals)), ', '));
        end
        placeholder = isnan(vals);                          % unknown lengths: generic values for the numerical cross-checks only
        gold = mod(0.6180339887*(1:numel(vals)), 1);
        vals(placeholder) = 0.1 + 0.4*gold(placeholder);
        if ~isempty(vars)
            real_vars = arrayfun(@(nm) sym(nm{1}, 'real'), names);
            rows = subs(rows, vars, real_vars);
        end
        dims = struct();
        for k = 1:numel(names), dims.(names{k}) = vals(k); end
        nm = matlab.lang.makeValidName(strtrim(edName.Value));
        mdh = mdh_table(nm, rows, dims);
        mdh.placeholder = placeholder;
    end

    function v = chooseVariants(hasDims)
        lens = ddLen.Value;  act = ddAct.Value;
        wantSym = ~strcmp(lens, 'Known (numeric values)');
        wantKnown = ~strcmp(lens, 'Unknown (symbolic)');
        if ~hasDims, wantSym = true; wantKnown = false; end   % without symbols the table is already numeric
        wantRigid = ~strcmp(act, 'With rotor inertia and friction');
        wantAct = ~strcmp(act, 'Rigid body only');
        v = {};
        if wantSym && wantRigid, v{end+1} = 'symbolic'; end
        if wantKnown && wantRigid, v{end+1} = 'known'; end
        if wantSym && wantAct, v{end+1} = 'actuation'; end
        if wantKnown && wantAct, v{end+1} = 'actuation_known'; end
    end

    function onComputeChanged()
        dynOnly = strcmp(ddCompute.Value, 'Dynamic model only (symbolic equations)');
        cbRec.Enable = ~dynOnly;  cbMin.Enable = ~dynOnly;
        btnRun.Text = ternary(dynOnly, 'RUN: derive the dynamic model', 'RUN: find and certify the base parameters');
    end

    function dynamicModelOnly()
        try
            mdh = buildMdh(false);                      % lengths may stay symbolic; not needed for the equations alone
            method = ddDeriv.Value;
            if mdh.n >= 7
                setStatus(sprintf('%d joints: symbolic derivation grows quickly; this can take several minutes.', mdh.n));
            else
                setStatus(sprintf('Deriving the dynamic model of %s by %s ...', mdh.name, method));
            end
            drawnow;
            t0 = tic;
            if strcmp(method, 'Newton-Euler')
                tau = mdh_rnea_symbolic(mdh);
            else
                tau = mdh_lagrange_torques(mdh);
            end
            dt = toc(t0);
            n = mdh.n;
            L = cell(n, 1);
            for j = 1:n
                kind = ternary(mdh.prismatic(j), 'force (N)', 'torque (N m)');
                L{j} = sprintf('tau_%d [%s] = %s', j, kind, dispExpr(tau(j)));
            end
            S.mdh = mdh;  S.out = [];  S.variants = {};  S.theta = struct();   % no base-parameter stage was run
            tblSum.Data = cell(0, 10);  ddVar.Items = {'-'};
            taParams.Value = L;
            lblCount.Text = sprintf('%d joint equations (%s), not regrouped into base parameters.', n, method);
            taLog.Value = {sprintf('Dynamic model only: %s, %d joints, %s derivation, %.1f s.', mdh.name, n, method, dt), ...
                'No monomial extraction, regrouping, or certification was run. Switch "What to compute" back to ' ...
                'find the base parameters, or re-enter this table there.'};
            setStatus(sprintf('Done in %.1f s. %d joint equations shown in the Results tab.', dt, n));
            tabs.SelectedTab = tabRes;
        catch err
            notify(err.message);
        end
    end

    function runPipeline()
        if strcmp(ddCompute.Value, 'Dynamic model only (symbolic equations)'), dynamicModelOnly();  return, end
        dlg = [];
        try
            needKnown = ~strcmp(ddLen.Value, 'Unknown (symbolic)');
            [~, vars0] = parseTable();
            if isempty(vars0), needKnown = false; end
            mdh = buildMdh(needKnown);
            variants = chooseVariants(~isempty(mdh.dim_names));
            if isempty(variants), notify('No variant selected.'); return, end
            opts = struct('results_root', edFolder.Value, 'variants', {variants}, 'parallel', cbPar.Value, ...
                'exact_reconstruction', cbRec.Value, 'exact_minimality', cbMin.Value, 'cache', true);
            setStatus(sprintf('Running %s (%d joints, %d prismatic), variants: %s ...', mdh.name, mdh.n, sum(mdh.prismatic), strjoin(variants, ', ')));
            drawnow;
            if strcmp(fig.Visible, 'on')
                dlg = uiprogressdlg(fig, 'Title', 'Running', 'Message', 'Deriving, regrouping and certifying (see the MATLAB command window) ...', 'Indeterminate', 'on');
            end
            t0 = tic;
            out = base_parameters_from_mdh(mdh, opts);
            if ~isempty(dlg) && isvalid(dlg), close(dlg); end
            S.mdh = mdh;  S.out = out;  S.variants = variants;
            showResults(out);
            msg = sprintf('Done in %.0f s. Results are in %s', toc(t0), out.folder);
            if any(mdh.placeholder)
                msg = [msg ' Unknown lengths were replaced by generic values only for the numerical cross-checks.'];
            end
            setStatus(msg);
            tabs.SelectedTab = tabRes;
        catch err
            if ~isempty(dlg) && isvalid(dlg), close(dlg); end
            notify(err.message);
        end
    end

    function showResults(out)
        T = out.table;  nT = numel(T);
        C = cell(nT, 10);
        for k = 1:nT
            C(k, :) = {T(k).variant, T(k).p, T(k).p_b, ternary(T(k).dependency_certified, 'certified', 'NOT certified'), ...
                T(k).reconstruction, T(k).minimality, T(k).e_max, T(k).r_geom, T(k).r_geom_rnea, round(T(k).time_s)};
        end
        tblSum.Data = C;
        ddVar.Items = {T.variant};  ddVar.Value = T(1).variant;
        showVariant();
        logf = fullfile(out.folder, 'run_log.txt');
        if isfile(logf), taLog.Value = splitlines(string(fileread(logf))); end
    end

    function th = loadTheta(variant)
        f = fullfile(S.out.folder, folder_of.(variant), 'base_regressor.mat');
        D = load(f, 'theta_base', 'A_matrix_base');
        th = D;
    end

    function showVariant()
        if isempty(S.out), return, end
        v = ddVar.Value;
        D = loadTheta(v);
        nb = numel(D.theta_base);
        L = cell(nb, 1);
        for k = 1:nb, L{k} = sprintf('%3d   %s', k, dispExpr(D.theta_base(k))); end
        taParams.Value = L;
        lblCount.Text = sprintf('%d base parameters', nb);
        S.theta.(v) = D;
    end

    function openFolder()
        if ~isempty(S.out) && isfolder(S.out.folder) && ispc, winopen(S.out.folder); end
    end

    function saveList()
        if isempty(S.out), return, end
        v = ddVar.Value;  f = fullfile(S.out.folder, folder_of.(v), 'theta_base_list.txt');
        setStatus(sprintf('The list is saved in %s', f));
        if strcmp(fig.Visible, 'on') && isfile(f) && ispc, winopen(fileparts(f)); end
    end

    function values = zeroActuationIfRigid(names, values, variant)
        % A "known"/"symbolic" variant excludes Ia, Fv, Fc from theta_base; zero them for the independent
        % Newton-Euler reference too, so the comparison matches what the chosen base-parameter set represents.
        if ~contains(variant, 'actuation')
            mask = ~cellfun('isempty', regexp(names, '^(Ia|Fv|Fc)\d+$', 'once'));
            values(mask) = 0;
        end
    end

    function [names, values] = standardValues(n, dims)
        M = tblStd.Data;  names = {};  values = [];
        for i = 1:n
            m = M(i, 1);  c = M(i, 2:4).';
            Ic = [M(i,5) M(i,8) M(i,9); M(i,8) M(i,6) M(i,10); M(i,9) M(i,10) M(i,7)];
            Io = Ic + m*((c.'*c)*eye(3) - c*c.');                          % parallel-axis theorem
            names = [names, {sprintf('m%d',i), sprintf('Lc%dx',i), sprintf('Lc%dy',i), sprintf('Lc%dz',i), ...
                sprintf('Ixx%d',i), sprintf('Iyy%d',i), sprintf('Izz%d',i), sprintf('Ixy%d',i), sprintf('Ixz%d',i), ...
                sprintf('Iyz%d',i), sprintf('Ia%d',i), sprintf('Fv%d',i), sprintf('Fc%d',i)}]; %#ok<AGROW>
            values = [values, m, c.', Io(1,1), Io(2,2), Io(3,3), Io(1,2), Io(1,3), Io(2,3), M(i,11), M(i,12), M(i,13)]; %#ok<AGROW>
        end
        if ~isempty(dims)
            names = [names, dims(:, 1).'];  values = [values, cell2mat(dims(:, 2)).'];
        end
    end

    function computeNominal()
        if isempty(S.out), notify('Run the pipeline first.'); return, end
        v = ddVar.Value;  D = S.theta.(v);
        n = S.mdh.n;  dims = tblDim.Data;
        if ~isempty(dims), dims = dims(cellfun(@(x) ~isempty(x) && isfinite(x), dims(:, 2)), :); end
        [names, values] = standardValues(n, dims);
        vars = symvar(D.theta_base);  vals = zeros(1, numel(vars));  missing = {};
        for j = 1:numel(vars)
            idx = find(strcmp(names, char(vars(j))), 1);
            if isempty(idx), missing{end+1} = char(vars(j)); else, vals(j) = values(idx); end %#ok<AGROW>
        end
        if ~isempty(missing)
            notify(['No value for: ' strjoin(missing, ', ') '. Give the link lengths (Robot tab) for variants with unknown lengths.']);
            return
        end
        th = double(subs(D.theta_base, vars, vals));
        C = cell(numel(th), 3);
        for k = 1:numel(th), C(k, :) = {k, dispExpr(D.theta_base(k)), th(k)}; end
        tblNom.Data = C;
        S.nominal = th;
        lblNom.Text = sprintf('%d numerical base parameters for the variant "%s".', numel(th), v);
    end

    function onLoadData()
        [f, p] = uigetfile({'*.mat;*.csv', 'Data files'}, 'Measured data');
        if isequal(f, 0), return, end
        try
            n = S.mdh.n;
            if endsWith(f, '.mat')
                Dd = load(fullfile(p, f));
                q = Dd.q;  qd = Dd.qd;  qdd = Dd.qdd;  tau = Dd.tau;
            else
                M = readmatrix(fullfile(p, f));
                if size(M, 2) ~= 4*n, error('The CSV file needs %d columns [q qd qdd tau], it has %d.', 4*n, size(M, 2)); end
                q = M(:, 1:n);  qd = M(:, n+1:2*n);  qdd = M(:, 2*n+1:3*n);  tau = M(:, 3*n+1:4*n);
            end
            identify(q, qd, qdd, tau);
        catch err
            notify(err.message);
        end
    end

    function res = identify(q, qd, qdd, tau)
        res = [];
        if isempty(S.out), notify('Run the pipeline first.'); return, end
        v = ddVar.Value;
        if ~ismember(v, {'known', 'actuation_known'})
            notify('Identification needs a variant with known lengths ("known" or "actuation_known").');  return
        end
        n = S.mdh.n;  N = size(q, 1);
        if any([size(qd, 1), size(qdd, 1), size(tau, 1)] ~= N) || any([size(q,2) size(qd,2) size(qdd,2) size(tau,2)] ~= n)
            notify(sprintf('q, qd, qdd and tau must all be N-by-%d.', n));  return
        end
        D = S.theta.(v);
        qs   = arrayfun(@(k) sym(sprintf('t%d',   k), 'real'), 1:n);
        qds  = arrayfun(@(k) sym(sprintf('t%dd',  k), 'real'), 1:n);
        qdds = arrayfun(@(k) sym(sprintf('t%ddd', k), 'real'), 1:n);
        Afun = matlabFunction(D.A_matrix_base, 'Vars', {qs, qds, qdds});
        nb = numel(D.theta_base);
        W = zeros(N*n, nb);  Y = zeros(N*n, 1);
        for k = 1:N
            W((k-1)*n + (1:n), :) = Afun(q(k, :), qd(k, :), qdd(k, :));
            Y((k-1)*n + (1:n)) = tau(k, :).';
        end
        th = W \ Y;
        res = struct('theta', th, 'rank', rank(W), 'cond', cond(W), 'rmse', sqrt(mean((W*th - Y).^2)), 'nb', nb);
        C = cell(nb, 3);
        for k = 1:nb, C(k, :) = {k, dispExpr(D.theta_base(k)), th(k)}; end
        tblId.Data = C;
        lblId.Text = sprintf('%d samples, %d base parameters; rank %d, condition number %.3g; torque RMSE %.3g.', N, nb, res.rank, res.cond, res.rmse);
        if res.rank < nb
            lblId2.Text = 'Warning: the data do not excite all base parameters (rank deficient). Use a richer trajectory.';
        else
            lblId2.Text = 'The data excite all base parameters (full column rank).';
        end
        S.ident = res;
    end

    function emax = rneaCheck()
        % Independent cross-check: the certified Lagrange regressor A_base*theta_base against the recursive
        % Newton-Euler algorithm (mdh_rnea_reference, Section 7.7 "Independent model"), at the same random
        % states and the same standard-parameter values. Both are exact identities in those values, so a large
        % disagreement signals an input error rather than round-off.
        emax = [];
        if isempty(S.mdh), notify('Run the pipeline first.'); return, end
        n = S.mdh.n;
        v = '';
        if ismember('known', S.variants), v = 'known';
        elseif ismember('actuation_known', S.variants), v = 'actuation_known'; end
        if isempty(v)
            notify('The Newton-Euler check needs a variant with known lengths ("known" or "actuation_known").');  return
        end
        if ~isfield(S.theta, v), S.theta.(v) = loadTheta(v); end
        D = S.theta.(v);
        dims = tblDim.Data;
        if ~isempty(dims), dims = dims(cellfun(@(x) ~isempty(x) && isfinite(x), dims(:, 2)), :); end
        [names, values] = standardValues(n, dims);
        if all(values(1:13:13*n) == 0)
            notify('Enter standard parameters (at least the link masses) in this tab first.');  return
        end
        values = zeroActuationIfRigid(names, values, v);
        vars = symvar(D.theta_base);  vals = zeros(1, numel(vars));  missing = {};
        for j = 1:numel(vars)
            idx = find(strcmp(names, char(vars(j))), 1);
            if isempty(idx), missing{end+1} = char(vars(j)); else, vals(j) = values(idx); end %#ok<AGROW>
        end
        if ~isempty(missing)
            notify(['No value for: ' strjoin(missing, ', ') '.']);  return
        end
        th = double(subs(D.theta_base, vars, vals));
        qs   = arrayfun(@(k) sym(sprintf('t%d',   k), 'real'), 1:n);
        qds  = arrayfun(@(k) sym(sprintf('t%dd',  k), 'real'), 1:n);
        qdds = arrayfun(@(k) sym(sprintf('t%ddd', k), 'real'), 1:n);
        Afun = matlabFunction(D.A_matrix_base, 'Vars', {qs, qds, qdds});
        fNE = mdh_rnea_reference(S.mdh);
        M = 5;  emax = 0;
        for s = 1:M
            q = 2*rand(1, n) - 1;  qd = 4*rand(1, n) - 2;  qdd = 6*rand(1, n) - 3;   % Table "Discovery ranges" of the paper
            tauL = (Afun(q, qd, qdd) * th).';
            tauN = fNE(q, qd, qdd, names, values).';
            emax = max(emax, max(abs(tauL - tauN)));
        end
        lblRnea.Text = sprintf(['Independent Newton-Euler check (variant "%s", %d random states): max |tau_{Lagrange} - ' ...
            'tau_{Newton-Euler}| = %.3g (N m or N). Agreement at this level confirms that the certified regressor ' ...
            'reproduces the full nonlinear dynamics, not only the entered standard parameters.'], v, M, emax);
        setStatus(sprintf('Newton-Euler check: max difference %.3g over %d random states.', emax, M));
    end

    function setTrajectory(duration_s, rate_hz, amplitude, noise_pct, outlier_pct)
        edDur.Value = duration_s;  edFs.Value = rate_hz;  edAmp.Value = amplitude;
        edNoise.Value = noise_pct;  edOut.Value = outlier_pct;
    end

    function res = generateAndIdentify()
        % Generates a quintic-ramped multisine excitation trajectory (same closed form as Eq. 6.4 of the paper),
        % simulates it with the independent Newton-Euler reference at the entered standard-parameter values, adds
        % the optional Gaussian noise and outliers of Section 7.9, and identifies the selected variant from it.
        res = [];
        if isempty(S.mdh), notify('Run the pipeline first.'); return, end
        n = S.mdh.n;
        if ~ismember(ddVar.Value, {'known', 'actuation_known'})
            notify('Identification needs a variant with known lengths ("known" or "actuation_known").');  return
        end
        dims = tblDim.Data;
        if ~isempty(dims), dims = dims(cellfun(@(x) ~isempty(x) && isfinite(x), dims(:, 2)), :); end
        [names, values] = standardValues(n, dims);
        if all(values(1:13:13*n) == 0)
            notify('Enter standard parameters (at least the link masses) in the Values tab first.');  return
        end
        values = zeroActuationIfRigid(names, values, ddVar.Value);
        Tdur = edDur.Value;  fs = edFs.Value;  Amp = edAmp.Value;  noisePct = edNoise.Value;  outPct = edOut.Value;
        Tr = min(5, Tdur/4);
        t = (0:1/fs:Tdur).';  N = numel(t);
        rho = min(t/Tr, 1);
        gam   = 10*rho.^3 - 15*rho.^4 + 6*rho.^5;
        dgam  = (30*rho.^2 - 60*rho.^3 + 30*rho.^4)/Tr;   dgam(t >= Tr)  = 0;
        ddgam = (60*rho - 180*rho.^2 + 120*rho.^3)/Tr^2;  ddgam(t >= Tr) = 0;
        q = zeros(N, n);  qd = q;  qdd = q;
        for i = 1:n
            fi = 0.05 + 0.01*(i - 1);  phi = (i - 1)*pi/5;  w = 2*pi*fi;  % distinct frequencies/phases, as in Eq. 6.4
            s = sin(w*t + phi);  c = cos(w*t + phi);
            q(:, i)   = Amp*gam.*s;
            qd(:, i)  = Amp*(dgam.*s + gam*w.*c);
            qdd(:, i) = Amp*(ddgam.*s + 2*dgam*w.*c - gam*w^2.*s);
        end
        fNE = mdh_rnea_reference(S.mdh);
        tau = zeros(N, n);
        for k = 1:N, tau(k, :) = fNE(q(k, :), qd(k, :), qdd(k, :), names, values).'; end
        if noisePct > 0
            sig = (noisePct/100)*std(tau, 0, 1);
            tau = tau + sig.*randn(N, n);
        end
        if outPct > 0
            sig = max((noisePct/100)*std(tau, 0, 1), eps);
            idxOut = randperm(N, min(round(outPct/100*N), N));
            signs = sign(randn(numel(idxOut), n));  signs(signs == 0) = 1;
            tau(idxOut, :) = tau(idxOut, :) + signs.*(15*sig);
        end
        res = identify(q, qd, qdd, tau);
        if ~isempty(res)
            setStatus(sprintf(['Generated %d samples (%.0f s at %.0f Hz) from the independent Newton-Euler reference; ' ...
                '%.0f%% noise, %.0f%% outliers. Identified below.'], N, Tdur, fs, noisePct, outPct));
        end
    end

    function saveTable(tb, defname)
        D = tb.Data;
        if isempty(D), return, end
        [f, p] = uiputfile('*.csv', 'Save as', defname);
        if isequal(f, 0), return, end
        T = cell2table(D, 'VariableNames', {'index', 'parameter', 'value'});
        writetable(T, fullfile(p, f));
        setStatus(sprintf('Saved %s', fullfile(p, f)));
    end

    function saveRobotFile()
        try
            [rows, vars] = parseTable();
        catch err
            notify(err.message);  return
        end
        nm = matlab.lang.makeValidName(strtrim(edName.Value));
        fname = fullfile(root, 'robots', ['mdh_' lower(nm) '.m']);
        if isfile(fname)
            if strcmp(fig.Visible, 'on')
                a = uiconfirm(fig, sprintf('%s exists. Overwrite?', fname), 'Save robot', 'Options', {'Overwrite', 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
                if strcmp(a, 'Cancel'), return, end
            else
                notify('The file exists; not overwritten.');  return
            end
        end
        D = tblMDH.Data;  n = size(D, 1);
        names = arrayfun(@char, vars, 'UniformOutput', false);
        L = {sprintf('function mdh = mdh_%s()', lower(nm)), ...
             sprintf('%%MDH_%s Modified DH table written by base_parameter_app.', upper(nm))};
        if ~isempty(names), L{end+1} = ['    syms ' strjoin(names, ' ') ' real']; end
        L{end+1} = '    rows = [';
        for i = 1:n
            L{end+1} = sprintf('        (%s), (%s), (%s), (%s), %d;', strtrim(D{i,2}), strtrim(D{i,3}), strtrim(D{i,4}), strtrim(D{i,5}), strcmp(D{i,1}, 'Prismatic')); %#ok<AGROW>
        end
        L{end+1} = '    ];';
        dd = tblDim.Data;  parts = {};
        for k = 1:numel(names)
            val = 1;
            if ~isempty(dd)
                j = find(strcmp(dd(:, 1), names{k}), 1);
                if ~isempty(j) && ~isempty(dd{j, 2}) && isfinite(dd{j, 2}), val = dd{j, 2}; end
            end
            parts{end+1} = sprintf('''%s'', %.10g', names{k}, val); %#ok<AGROW>
        end
        if isempty(parts), L{end+1} = '    dims = struct();'; else, L{end+1} = ['    dims = struct(' strjoin(parts, ', ') ');']; end
        L{end+1} = sprintf('    mdh = mdh_table(''%s'', rows, dims);', nm);
        L{end+1} = 'end';
        fid = fopen(fname, 'w');  fprintf(fid, '%s\n', L{:});  fclose(fid);
        ddPreset.Items = presetList();
        setStatus(sprintf('Saved %s. It can be used in scripts as mdh = mdh_%s();', fname, lower(nm)));
    end

    function out = ternary(c, a, b)
        if c, out = a; else, out = b; end
    end
end
