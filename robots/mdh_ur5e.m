function mdh = mdh_ur5e()
%MDH_UR5E Modified DH table of the Universal Robots UR5e (manufacturer zero configuration
%   shifted into the MDH frames; see the paper, Table 4). Columns: [alpha_{i-1} a_{i-1} d_i theta0_i].
    syms L1 L2 L3 L4 L5 L6 LL1 LL2 real
    rows = [ 0      0    L1    -pi;
             pi/2   0    LL1    pi/2;
             0      L2  -LL2    0;
            -pi     L3   L4    -pi/2;
            -pi/2   0    L5     0;
             pi/2   0    L6     pi/2 ];
    dims = struct('L1', 0.1625, 'L2', 0.425, 'L3', 0.3922, 'L4', -0.1333, ...
                  'L5', -0.0997, 'L6', -0.0996, 'LL1', 0, 'LL2', 0);
    mdh = mdh_table('UR5e', rows, dims);
end
