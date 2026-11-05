function mdh = mdh_scara_up()
%MDH_SCARA_UP The SCARA of mdh_scara with the prismatic axis drawn upward (alpha_2 = 0 instead of pi).
%   The robot is the same; only the sign of the prismatic coordinate q3 and of the tool joint q4 changes.
%   The base-parameter count must therefore be unchanged, which this table checks.
%   Columns: [alpha_{i-1} a_{i-1} d_i theta0_i type]; type 1 = prismatic.
    syms L1 L2 L3 L4 real
    rows = [ 0     0    L1   0   0;
             0     L2   0    0   0;
             0     L3   L4   0   1;
             0     0    0    0   0 ];
    dims = struct('L1', 0.40, 'L2', 0.325, 'L3', 0.275, 'L4', 0.10);
    mdh = mdh_table('SCARA_up', rows, dims);
end
