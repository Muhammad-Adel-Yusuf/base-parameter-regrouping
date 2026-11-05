function P = extract_nominal_parameters(robot_name, bodies, mdh, S, delta)
%EXTRACT_NOMINAL_PARAMETERS Link inertial data of a loadrobot model in MDH frames.
%
%   P = extract_nominal_parameters(robot_name, bodies, mdh, S, delta)
%
%   bodies : names of the n moving bodies (link 1..n)
%   mdh    : struct with fields alpha (1xn), a (1xn), d (1xn), theta0 (1xn)
%            describing rows 1..n of the MDH table (theta_i = q_i + theta0_i)
%   S,delta: joint mapping q_mdh = S.*q_m + delta
%
%   Returns P.m (1xn), P.Lc (nx3) centre of mass in MDH frame {i},
%   P.Ic (nx6) centroidal inertia [Ixx Iyy Izz Ixy Ixz Iyz] in the axes of
%   MDH frame {i} (tensor elements), and P.attach_error, the largest change
%   of the body-to-MDH transform over random configurations (must be ~0).

    r = panda_robot();  n = numel(bodies);   % robot_name: loadrobot id, for the record
    rng(11);  Q = (rand(6, n) - 0.5)*3;
    P.m = zeros(1,n); P.Lc = zeros(n,3); P.Ic = zeros(n,6); P.attach_error = 0;
    for i = 1:n
        b = getBody(r, bodies{i});
        for k = 1:size(Q,1)
            qm = Q(k,:);  q = S.*qm + delta;
            Tm = mdh_frame(mdh, q, i);
            X = Tm \ getTransform(r, qm, bodies{i});
            if k == 1, X0 = X; else, P.attach_error = max(P.attach_error, norm(X - X0, 'fro')); end
        end
        R = X0(1:3,1:3);  p = X0(1:3,4);
        m = b.Mass;  c = b.CenterOfMass(:);
        J = b.Inertia;   % [Ixx Iyy Izz Iyz Ixz Ixy] about the body origin, body axes
        Io = [J(1) J(6) J(5); J(6) J(2) J(4); J(5) J(4) J(3)];
        Ic_body = Io - m*((c.'*c)*eye(3) - c*c.');     % parallel-axis shift to the CoM
        Ic = R*Ic_body*R.';                              % rotate into MDH axes
        P.m(i) = m;
        P.Lc(i,:) = (R*c + p).';
        P.Ic(i,:) = [Ic(1,1) Ic(2,2) Ic(3,3) Ic(1,2) Ic(1,3) Ic(2,3)];
    end
    % remove floating-point residue of the frame transformations
    P.Lc(abs(P.Lc) < 1e-9) = 0;
    P.Ic(abs(P.Ic) < 1e-9) = 0;
end

function T = mdh_frame(mdh, q, i)
    T = eye(4);
    for j = 1:i
        al = mdh.alpha(j); a = mdh.a(j); d = mdh.d(j); th = q(j) + mdh.theta0(j);
        T = T*[cos(th) -sin(th) 0 a;
               sin(th)*cos(al) cos(th)*cos(al) -sin(al) -d*sin(al);
               sin(th)*sin(al) cos(th)*sin(al)  cos(al)  d*cos(al);
               0 0 0 1];
    end
end
