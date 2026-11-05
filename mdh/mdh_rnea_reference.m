function f = mdh_rnea_reference(mdh)
%MDH_RNEA_REFERENCE Independent numerical inverse dynamics (recursive Newton-Euler).
%
%   f = mdh_rnea_reference(mdh)
%   tau = f(q, qd, qdd, names, values)
%
%   Returns a function handle that evaluates the recursive Newton-Euler
%   algorithm for the MDH table with origin-referenced inertia parameters.
%   It does not use the symbolic torque equations and serves as the independent
%   reference of the verification stage (Stage 5) of the regrouping. The
%   physical parameters are given by name (m1, Lc1x, Ixx1, ..., link dimensions,
%   Ia1, Fv1, Fc1); a parameter that is not listed is zero.

    n = mdh.n;
    T = [mdh.alpha(:) mdh.a(:) mdh.d(:) mdh.theta0(:)];
    dims = mdh.dim_symbols;
    if isempty(dims)
        tab = @(x) double(T);
    else
        tab = matlabFunction(T, 'Vars', {dims});
    end
    dim_names = mdh.dim_names;
    pr = mdh.prismatic(:);
    f = @(q, qd, qdd, names, values) rnea(tab, dim_names, pr, n, q, qd, qdd, names, values);
end

function tau = rnea(tab, dim_names, pr, n, q, qd, qdd, names, values)
    v = containers.Map(cellstr(string(names)), num2cell(double(values)));
    get = @(s) lookup(v, s);
    dv = zeros(1, numel(dim_names));
    for k = 1:numel(dim_names), dv(k) = get(dim_names{k}); end
    TT = tab(dv);
    alpha = TT(:,1);  a = TT(:,2);  d = TT(:,3);  th0 = TT(:,4);
    R = cell(1, n+1);  p = cell(1, n+1);
    for i = 1:n
        if pr(i), th = th0(i);  di = d(i) + q(i);          % prismatic: q adds to the offset
        else,     th = q(i) + th0(i);  di = d(i);           % revolute: q adds to the angle
        end
        ca = cos(alpha(i)); sa = sin(alpha(i)); ct = cos(th); st = sin(th);
        R{i} = [ct -st 0; st*ca ct*ca -sa; st*sa ct*sa ca];
        p{i} = [a(i); -di*sa; di*ca];
    end
    R{n+1} = eye(3);  p{n+1} = zeros(3,1);

    z = [0;0;1];
    w = zeros(3,1);  dw = zeros(3,1);  acc = [0;0;9.81];
    F = cell(1,n);  N = cell(1,n);
    for i = 1:n
        Rt = R{i}.';
        acc = Rt*(acc + cross(dw, p{i}) + cross(w, cross(w, p{i})));
        dw  = Rt*dw;
        w   = Rt*w;
        if pr(i)
            acc = acc + qdd(i)*z + 2*cross(w, qd(i)*z);     % sliding joint: Coriolis term
        else
            dw  = dw + cross(w, qd(i)*z) + qdd(i)*z;
            w   = w + qd(i)*z;
        end
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
        if pr(i), tau(i) = f(3); else, tau(i) = nn(3); end     % force along z_i or moment about z_i
    end
    for i = 1:n                                        % joint-side actuation terms, if present
        tau(i) = tau(i) + get(sprintf('Ia%d', i))*qdd(i) + get(sprintf('Fv%d', i))*qd(i) ...
                        + get(sprintf('Fc%d', i))*sign(qd(i));
    end
end

function x = lookup(v, s)
    if isKey(v, s), x = v(s); else, x = 0; end
end
