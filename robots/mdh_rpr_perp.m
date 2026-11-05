function mdh = mdh_rpr_perp()
%MDH_RPR_PERP Revolute-prismatic-revolute test robot with perpendicular axes (see mdh_rpr_twist).
%   Columns: [alpha_{i-1} a_{i-1} d_i theta0_i type].
    syms L1 L2 L3 L4 L5 real
    rows = [ 0      0    L1   0   0;
             pi/2   L2   L3   0   1;
             pi/2   L4   L5   0   0 ];
    dims = struct('L1', 0.30, 'L2', 0.20, 'L3', 0.15, 'L4', 0.25, 'L5', 0.10);
    mdh = mdh_table('RPR_perp', rows, dims);
end
