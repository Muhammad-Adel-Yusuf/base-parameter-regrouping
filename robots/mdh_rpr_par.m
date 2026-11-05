function mdh = mdh_rpr_par()
%MDH_RPR_PAR Revolute-prismatic-revolute test robot with parallel axes (see mdh_rpr_twist).
%   Columns: [alpha_{i-1} a_{i-1} d_i theta0_i type].
    syms L1 L2 L3 L4 L5 real
    rows = [ 0   0    L1   0   0;
             0   L2   L3   0   1;
             0   L4   L5   0   0 ];
    dims = struct('L1', 0.30, 'L2', 0.20, 'L3', 0.15, 'L4', 0.25, 'L5', 0.10);
    mdh = mdh_table('RPR_par', rows, dims);
end
