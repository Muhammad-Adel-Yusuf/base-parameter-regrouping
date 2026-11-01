function [robot, P] = panda_robot()
%PANDA_ROBOT Franka Emika Panda arm with published inertial parameters.
%
%   [robot, P] = panda_robot()
%
% The kinematics and visual meshes are those of loadrobot('frankaEmikaPanda')
% (Franka Emika description, Apache-2.0). That description carries no mass
% properties, so the masses, centres of mass and centroidal inertia tensors
% of links 1..7 are set to the feasible parameters identified by
%   C. Gaz, M. Cognetti, A. Oliva, P. Robuffo Giordano, A. De Luca,
%   "Dynamic identification of the Franka Emika Panda robot with retrieval
%   of feasible parameters using penalty-based optimization",
%   IEEE RA-L 4(4):4147-4154, 2019, supplementary material, Tables VIII-IX.
% The hand and fingers are removed (arm without end effector). The body
% frames of panda_link1..7 coincide with the modified DH (Craig) frames of
% A_Full_Dynamics_Panda_parallel_axis.m, so the values below are expressed
% in the MDH frames directly. c1z does not enter the joint torques (it is
% reported as "don't care" by Gaz et al.); the value of the alternative
% identification in the same table is used.
%
%   P.m (1x7), P.Lc (7x3) centre of mass in frame {i} (m),
%   P.Ic (7x6) centroidal inertia [Ixx Iyy Izz Ixy Ixz Iyz] (kg m^2).

    P.m  = [4.970684 0.646926 3.228604 3.587895 1.225946 1.666555 0.735522];
    P.Lc = [ 3.875e-03   2.081e-03  -0.1750
            -3.141e-03  -2.872e-02   3.495e-03
             2.7518e-02  3.9252e-02 -6.6502e-02
            -5.317e-02   1.04419e-01 2.7454e-02
            -1.1953e-02  4.1065e-02 -3.8437e-02
             6.0149e-02 -1.4117e-02 -1.0517e-02
             1.0517e-02 -4.252e-03   6.1597e-02];
    %       Ixx         Iyy         Izz         Ixy         Ixz         Iyz
    P.Ic = [7.0337e-01  7.0661e-01  9.1170e-03 -1.3900e-04  6.7720e-03  1.9169e-02
            7.9620e-03  2.8110e-02  2.5995e-02 -3.9250e-03  1.0254e-02  7.0400e-04
            3.7242e-02  3.6155e-02  1.0830e-02 -4.7610e-03 -1.1396e-02 -1.2805e-02
            2.5853e-02  1.9552e-02  2.8323e-02  7.7960e-03 -1.3320e-03  8.6410e-03
            3.5549e-02  2.9474e-02  8.6270e-03 -2.1170e-03 -4.0370e-03  2.2900e-04
            1.9640e-03  4.3540e-03  5.4330e-03  1.0900e-04 -1.1580e-03  3.4100e-04
            1.2516e-02  1.0027e-02  4.8150e-03 -4.2800e-04 -1.1960e-03 -7.4100e-04];

    robot = loadrobot('frankaEmikaPanda', 'DataFormat', 'row');
    removeBody(robot, 'panda_hand');            % also removes both fingers
    robot.Gravity = [0 0 -9.81];
    for i = 1:7
        b = getBody(robot, sprintf('panda_link%d', i));
        c  = P.Lc(i,:).';
        Ic = [P.Ic(i,1) P.Ic(i,4) P.Ic(i,5); P.Ic(i,4) P.Ic(i,2) P.Ic(i,6); P.Ic(i,5) P.Ic(i,6) P.Ic(i,3)];
        if min(eig(Ic)) <= 0, error('Link %d: centroidal inertia is not positive definite.', i); end
        Io = Ic + P.m(i)*((c.'*c)*eye(3) - c*c.');            % about the body origin
        b.Mass = P.m(i);
        b.CenterOfMass = c.';
        b.Inertia = [Io(1,1) Io(2,2) Io(3,3) Io(2,3) Io(1,3) Io(1,2)];
    end
    b = getBody(robot, 'panda_link8');           % flange frame, massless
    b.Mass = 0;  b.Inertia = zeros(1,6);

    % Visual meshes only (no effect on the dynamics): the description shipped
    % with MATLAB uses low-polygon meshes, so the detailed visual meshes of
    % franka_description (Apache-2.0; simulink/meshes_franka_visual/linkK.stl,
    % converted from the COLLADA files by simulink/dae_to_stl.m) replace them.
    mesh_dir = fullfile(fileparts(mfilename('fullpath')), 'simulink', 'meshes_franka_visual');
    for i = 0:7
        f = fullfile(mesh_dir, sprintf('link%d.stl', i));
        if ~isfile(f), continue; end
        if i == 0, b = robot.Base; else, b = getBody(robot, sprintf('panda_link%d', i)); end
        clearVisual(b);
        addVisual(b, 'Mesh', f);
    end
end
