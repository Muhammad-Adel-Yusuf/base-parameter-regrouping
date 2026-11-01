function tau = panda_rnea_numeric(q, qd, qdd, names, values)
%PANDA_RNEA_NUMERIC Numerical inverse dynamics of the Panda (origin-referenced inertia).
%
%   tau = panda_rnea_numeric(q, qd, qdd, names, values)
%
%   Recursive Newton-Euler algorithm with the same modified DH table and
%   parameterization as A_Full_Dynamics_Panda_parallel_axis.m. The physical
%   parameters are given by name (m1, Lc1x, Ixx1, ..., L1..L6); a parameter
%   that is not listed is zero. Used as an independent reference for the
%   verification stage (Stage 5) of the regrouping.

    v = containers.Map(cellstr(string(names)), num2cell(double(values)));
    get = @(s) lookup(v, s);
    L = arrayfun(@(k) get(sprintf('L%d', k)), 1:6);
    alpha = [0 -pi/2 pi/2 pi/2 -pi/2 pi/2 pi/2];
    a     = [0 0 0 L(3) -L(4) 0 L(6)];
    d     = [L(1) 0 L(2) 0 L(5) 0 0];
    n = 7;
    R = cell(1, n+1);  p = cell(1, n+1);
    for i = 1:n
        ca = cos(alpha(i)); sa = sin(alpha(i)); ct = cos(q(i)); st = sin(q(i));
        R{i} = [ct -st 0; st*ca ct*ca -sa; st*sa ct*sa ca];
        p{i} = [a(i); -d(i)*sa; d(i)*ca];
    end
    R{n+1} = eye(3);  p{n+1} = zeros(3,1);

    z = [0;0;1];
    w = zeros(3,1);  dw = zeros(3,1);  acc = [0;0;9.81];
    F = cell(1,n);  N = cell(1,n);
    for i = 1:n
        Rt = R{i}.';
        acc = Rt*(acc + cross(dw, p{i}) + cross(w, cross(w, p{i})));
        dw  = Rt*dw + cross(Rt*w, qd(i)*z) + qdd(i)*z;
        w   = Rt*w + qd(i)*z;
        m  = get(sprintf('m%d', i));
        MS = m*[get(sprintf('Lc%dx', i)); get(sprintf('Lc%dy', i)); get(sprintf('Lc%dz', i))];
        Ixx = get(sprintf('Ixx%d',i)); Iyy = get(sprintf('Iyy%d',i)); Izz = get(sprintf('Izz%d',i));
        Ixy = get(sprintf('Ixy%d',i)); Ixz = get(sprintf('Ixz%d',i)); Iyz = get(sprintf('Iyz%d',i));
        I = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];
        F{i} = m*acc + cross(dw, MS) + cross(w, cross(w, MS));
        N{i} = I*dw + cross(w, I*w) + cross(MS, acc);
    end
    f = zeros(3,1);  nn = zeros(3,1);  tau = zeros(n,1);
    for i = n:-1:1
        nn = N{i} + R{i+1}*nn + cross(p{i+1}, R{i+1}*f);
        f  = F{i} + R{i+1}*f;
        tau(i) = nn(3);
    end
    % joint-side actuation terms, if present
    for i = 1:n
        tau(i) = tau(i) + get(sprintf('Ia%d', i))*qdd(i) + get(sprintf('Fv%d', i))*qd(i) ...
                        + get(sprintf('Fc%d', i))*sign(qd(i));
    end
end

function x = lookup(v, s)
    if isKey(v, s), x = v(s); else, x = 0; end
end
