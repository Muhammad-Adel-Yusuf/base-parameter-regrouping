function [tau, q, qd, qdd] = mdh_lagrange_torques(mdh, g_val)
%MDH_LAGRANGE_TORQUES Symbolic Euler-Lagrange inverse dynamics from an MDH table.
%
%   [tau, q, qd, qdd] = mdh_lagrange_torques(mdh, g_val)
%
%   mdh    struct from mdh_table.m
%   g_val  gravitational acceleration (default sym(981)/100)
%
%   The energies use origin-referenced (parallel-axis) inertia tensors,
%       K_i = 1/2 m_i v_i'v_i + m_i v_i'(w_i x c_i) + 1/2 w_i' I_oi w_i
%       V_i = -m_i g0'(R_i c_i + o_i),
%   with the origin velocity v_i and angular velocity w_i of frame {i}
%   propagated in {i}; the torques are linear in the 10 standard parameters of
%   every link. For a prismatic joint the entry tau(j) is the generalized force
%   along z_j (a force, not a torque) and q(j) is the displacement. The kinematic symbols are t<k> (position), t<k>d (velocity)
%   and t<k>dd (acceleration), k = 1..n, as in the regrouping code.
%   The link dimensions stay symbolic. The expressions are returned
%   unsimplified (nested): the monomials are extracted afterwards.
%
%   Inertial symbols: m<i>, Lc<i>x/y/z (centre of mass, frame {i}), and the
%   origin-referenced tensor elements Ixx<i>, Iyy<i>, Izz<i>, Ixy<i>, Ixz<i>, Iyz<i>.

    if nargin < 2, g_val = sym(981)/100; end
    n = mdh.n;
    q   = arrayfun(@(k) sym(sprintf('t%d',   k), 'real'), 1:n);
    qd  = arrayfun(@(k) sym(sprintf('t%dd',  k), 'real'), 1:n);
    qdd = arrayfun(@(k) sym(sprintf('t%ddd', k), 'real'), 1:n);
    m = sym('m', [1 n], 'real');

    w = sym(zeros(3,1));  v = sym(zeros(3,1));  T0 = sym(eye(4));
    KE = sym(0);  PE = sym(0);
    for i = 1:n
        al = mdh.alpha(i);  a = mdh.a(i);
        if mdh.prismatic(i)
            d = q(i) + mdh.d(i);  th = mdh.theta0(i);          % prismatic: q_i adds to the offset
        else
            d = mdh.d(i);  th = q(i) + mdh.theta0(i);          % revolute: q_i adds to the angle
        end
        Ri = [cos(th)          -sin(th)          0;
              sin(th)*cos(al)   cos(th)*cos(al) -sin(al);
              sin(th)*sin(al)   cos(th)*sin(al)  cos(al)];
        pim1_i = [a; -d*sin(al); d*cos(al)];
        Lc = [sym(sprintf('Lc%dx',i),'real'); sym(sprintf('Lc%dy',i),'real'); sym(sprintf('Lc%dz',i),'real')];
        Ixx = sym(sprintf('Ixx%d',i),'real'); Iyy = sym(sprintf('Iyy%d',i),'real'); Izz = sym(sprintf('Izz%d',i),'real');
        Ixy = sym(sprintf('Ixy%d',i),'real'); Ixz = sym(sprintf('Ixz%d',i),'real'); Iyz = sym(sprintf('Iyz%d',i),'real');
        I = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];
        v = Ri.'*(v + cross(w, pim1_i));
        w = Ri.'*w;
        if mdh.prismatic(i)
            v = v + [0;0;qd(i)];                               % sliding along z_i
        else
            w = w + [0;0;qd(i)];
        end
        T0 = T0*[Ri pim1_i; 0 0 0 1];
        KE = KE + m(i)/2*(v.'*v) + m(i)*(v.'*cross(w, Lc)) + (w.'*I*w)/2;
        PE = PE + m(i)*[0 0 g_val]*(T0(1:3,1:3)*Lc + T0(1:3,4));
    end
    Lag = KE - PE;
    tau = sym(zeros(n,1));
    for j = 1:n
        P = diff(Lag, qd(j));
        tau(j) = sum(arrayfun(@(k) diff(P, q(k))*qd(k) + diff(P, qd(k))*qdd(k), 1:n)) - diff(Lag, q(j));
    end
end
