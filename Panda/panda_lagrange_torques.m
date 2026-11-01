function tau = panda_lagrange_torques(q, qd, qdd, L, g_val)
%PANDA_LAGRANGE_TORQUES Symbolic Euler-Lagrange inverse dynamics of the Panda.
%
%   tau = panda_lagrange_torques(q, qd, qdd, L, g_val)
%
%   q, qd, qdd : 1x7 symbolic joint coordinates, velocities, accelerations
%   L          : 1x6 link dimensions (symbols or numbers), d1 = L(1), d3 = L(2),
%                a3 = L(3), a4 = -L(4), d5 = L(5), a6 = L(6)
%   g_val      : gravitational acceleration
%
%   The energies use origin-referenced inertia tensors (parallel-axis form), as
%   in the UR5e and IRB 120 scripts,
%       K_i = 1/2 m_i v_i'v_i + m_i v_i'(w_i x c_i) + 1/2 w_i' I_oi w_i
%       V_i = -m_i g0'(R_i c_i + o_i),
%   with the origin velocity v_i and angular velocity w_i of frame {i} propagated
%   in {i}. The torques follow from tau_j = d/dt(dL/dqd_j) - dL/dq_j, where the
%   total time derivative is the chain rule
%       d/dt f(q,qd) = sum_k [ df/dq_k qd_k + df/dqd_k qdd_k ].
%   The expressions are returned unsimplified (nested), which keeps the
%   derivation to seconds; the monomials are extracted afterwards.

    nL = 7;
    DH = [ 0      0     q(1)  L(1);
          -pi/2   0     q(2)  0;
           pi/2   0     q(3)  L(2);
           pi/2   L(3)  q(4)  0;
          -pi/2  -L(4)  q(5)  L(5);
           pi/2   0     q(6)  0;
           pi/2   L(6)  q(7)  0 ];
    m = sym('m', [1 nL], 'real');
    w = sym(zeros(3,1));  v = sym(zeros(3,1));  T0 = sym(eye(4));
    KE = sym(0);  PE = sym(0);
    for i = 1:nL
        al = DH(i,1); a = DH(i,2); th = DH(i,3); d = DH(i,4);
        Ri = [cos(th)          -sin(th)          0;
              sin(th)*cos(al)   cos(th)*cos(al) -sin(al);
              sin(th)*sin(al)   cos(th)*sin(al)  cos(al)];
        pim1_i = [a; -d*sin(al); d*cos(al)];
        Lc = [sym(sprintf('Lc%dx',i),'real'); sym(sprintf('Lc%dy',i),'real'); sym(sprintf('Lc%dz',i),'real')];
        Ixx = sym(sprintf('Ixx%d',i),'real'); Iyy = sym(sprintf('Iyy%d',i),'real'); Izz = sym(sprintf('Izz%d',i),'real');
        Ixy = sym(sprintf('Ixy%d',i),'real'); Ixz = sym(sprintf('Ixz%d',i),'real'); Iyz = sym(sprintf('Iyz%d',i),'real');
        I = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];            % origin-referenced tensor
        v = Ri.'*(v + cross(w, pim1_i));
        w = Ri.'*w + [0;0;qd(i)];
        T0 = T0*[Ri pim1_i; 0 0 0 1];
        KE = KE + m(i)/2*(v.'*v) + m(i)*(v.'*cross(w, Lc)) + (w.'*I*w)/2;
        PE = PE + m(i)*[0 0 g_val]*(T0(1:3,1:3)*Lc + T0(1:3,4));
    end
    Lag = KE - PE;
    tau = sym(zeros(nL,1));
    for j = 1:nL
        P = diff(Lag, qd(j));
        tau(j) = sum(arrayfun(@(k) diff(P, q(k))*qd(k) + diff(P, qd(k))*qdd(k), 1:nL)) - diff(Lag, q(j));
    end
end
