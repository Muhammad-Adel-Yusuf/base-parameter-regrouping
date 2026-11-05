function mdh = mdh_scara()
%MDH_SCARA Modified DH table of a SCARA robot (RRPR): two horizontal revolute joints, a vertical
%   prismatic joint and a revolute tool joint whose axis is collinear with the prismatic one.
%   Columns: [alpha_{i-1} a_{i-1} d_i theta0_i type]; type 1 = prismatic (q3 is added to d3 = L4).
%   z3 points downwards (alpha_2 = pi), the usual SCARA convention.
    syms L1 L2 L3 L4 real
    rows = [ 0     0    L1   0   0;
             0     L2   0    0   0;
             pi    L3   L4   0   1;
             0     0    0    0   0 ];
    dims = struct('L1', 0.40, 'L2', 0.325, 'L3', 0.275, 'L4', 0.10);
    mdh = mdh_table('SCARA', rows, dims);
end
