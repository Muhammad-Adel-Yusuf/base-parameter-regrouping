function [tau, q, qd, qdd] = mdh_rnea_symbolic(mdh, g_val)
%MDH_RNEA_SYMBOLIC Symbolic recursive Newton-Euler inverse dynamics from an MDH table.
%
%   [tau, q, qd, qdd] = mdh_rnea_symbolic(mdh, g_val)
%
%   mdh    struct from mdh_table.m
%   g_val  gravitational acceleration (default sym(981)/100)
%
%   Same recursive Newton-Euler algorithm as the independent numerical reference
%   (mdh_rnea_reference.m), evaluated symbolically instead of substituting numerical
%   parameter values: the base acceleration is set to g (upward) to account for gravity,
%   angular and linear quantities are propagated outward link by link, and forces and
%   moments are propagated back from the tip. For a prismatic joint tau(j) is the
%   generalized force along z_j and q(j) is the displacement.
%
%   Kinematic symbols t<k>, t<k>d, t<k>dd and inertial symbols m<i>, Lc<i>x/y/z,
%   Ixx<i>, Iyy<i>, Izz<i>, Ixy<i>, Ixz<i>, Iyz<i> are named exactly as in
%   mdh_lagrange_torques.m, so the two derivations can be compared term by term or
%   substituted with the same numerical values. Rotor inertia and friction are not
%   included here (added later, as in the Euler-Lagrange route). The link dimensions
%   stay symbolic; the expressions are returned unsimplified (nested).

    if nargin < 2, g_val = sym(981)/100; end
    n = mdh.n;
    q   = arrayfun(@(k) sym(sprintf('t%d',   k), 'real'), 1:n);
    qd  = arrayfun(@(k) sym(sprintf('t%dd',  k), 'real'), 1:n);
    qdd = arrayfun(@(k) sym(sprintf('t%ddd', k), 'real'), 1:n);

    R = cell(1, n + 1);  p = cell(1, n + 1);
    for i = 1:n
        al = mdh.alpha(i);  a = mdh.a(i);
        if mdh.prismatic(i)
            d = q(i) + mdh.d(i);  th = mdh.theta0(i);          % prismatic: q_i adds to the offset
        else
            d = mdh.d(i);  th = q(i) + mdh.theta0(i);          % revolute: q_i adds to the angle
        end
        R{i} = [cos(th)          -sin(th)          0;
                sin(th)*cos(al)   cos(th)*cos(al) -sin(al);
                sin(th)*sin(al)   cos(th)*sin(al)  cos(al)];
        p{i} = [a; -d*sin(al); d*cos(al)];
    end
    R{n+1} = sym(eye(3));  p{n+1} = sym(zeros(3,1));

    z = sym([0; 0; 1]);
    w = sym(zeros(3,1));  dw = sym(zeros(3,1));  acc = sym([0; 0; g_val]);
    F = cell(1, n);  N = cell(1, n);
    for i = 1:n
        Rt = R{i}.';
        acc = Rt*(acc + cross(dw, p{i}) + cross(w, cross(w, p{i})));
        dw  = Rt*dw;
        w   = Rt*w;
        if mdh.prismatic(i)
            acc = acc + qdd(i)*z + 2*cross(w, qd(i)*z);        % sliding joint: Coriolis term
        else
            dw  = dw + cross(w, qd(i)*z) + qdd(i)*z;
            w   = w + qd(i)*z;
        end
        m  = sym(sprintf('m%d', i), 'real');
        Lc = [sym(sprintf('Lc%dx',i),'real'); sym(sprintf('Lc%dy',i),'real'); sym(sprintf('Lc%dz',i),'real')];
        MS = m*Lc;
        Ixx = sym(sprintf('Ixx%d',i),'real'); Iyy = sym(sprintf('Iyy%d',i),'real'); Izz = sym(sprintf('Izz%d',i),'real');
        Ixy = sym(sprintf('Ixy%d',i),'real'); Ixz = sym(sprintf('Ixz%d',i),'real'); Iyz = sym(sprintf('Iyz%d',i),'real');
        I = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];
        F{i} = m*acc + cross(dw, MS) + cross(w, cross(w, MS));
        N{i} = I*dw + cross(w, I*w) + cross(MS, acc);
    end
    f = sym(zeros(3,1));  nn = sym(zeros(3,1));  tau = sym(zeros(n,1));
    for i = n:-1:1
        nn = N{i} + R{i+1}*nn + cross(p{i+1}, R{i+1}*f);
        f  = F{i} + R{i+1}*f;
        if mdh.prismatic(i), tau(i) = f(3); else, tau(i) = nn(3); end
    end
end
