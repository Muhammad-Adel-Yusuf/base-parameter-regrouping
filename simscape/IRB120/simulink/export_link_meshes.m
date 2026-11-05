function files = export_link_meshes(robot, folder, prefix)
%EXPORT_LINK_MESHES Write the visual mesh of every body of a rigidBodyTree
% to an STL file expressed in that body's frame.
%
%   files = export_link_meshes(robot, folder, prefix)
%
% The meshes of the robot descriptions shipped with MATLAB (loadrobot) are
% stored inside the rigidBodyTree object. They are recovered from the
% patches drawn by show() at the home configuration and transformed back
% into the body frames. files is a struct array with fields body and file.
% The meshes originate from the ROS-Industrial / Universal Robots robot
% descriptions (BSD / Apache-2.0 licences, see the .rights files in
% toolbox/robotics/robotmodels/robotModels/roboturdf).

    if ~isfolder(folder), mkdir(folder); end
    q0 = homeConfiguration(robot);
    fig = figure('Visible', 'off');  cleanup = onCleanup(@() close(fig));
    ax = show(robot, q0, 'Frames', 'off', 'Visuals', 'on');
    patches = findobj(ax, 'Type', 'patch');
    files = struct('body', {}, 'file', {});
    names = [{robot.BaseName}, cellfun(@(b) b.Name, robot.Bodies, 'UniformOutput', false)];
    for k = 1:numel(names)
        p = findobj(patches, 'DisplayName', [names{k} '_mesh']);
        if isempty(p), continue; end
        T = getTransform(robot, q0, names{k});
        V = double([p(1).Vertices, ones(size(p(1).Vertices, 1), 1)]) / T.';   % world -> body frame
        f = fullfile(folder, sprintf('%s_%s.stl', prefix, names{k}));
        stlwrite(triangulation(double(p(1).Faces), V(:, 1:3)), f);
        files(end+1) = struct('body', names{k}, 'file', f); %#ok<AGROW>
    end
end
