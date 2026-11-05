function mdh = mdh_panda()
%MDH_PANDA Modified DH table of the Franka Emika Panda (arm without hand).
%   Columns: [alpha_{i-1} a_{i-1} d_i theta0_i].
    syms L1 L2 L3 L4 L5 L6 real
    rows = [ 0      0    L1   0;
            -pi/2   0    0    0;
             pi/2   0    L2   0;
             pi/2   L3   0    0;
            -pi/2  -L4   L5   0;
             pi/2   0    0    0;
             pi/2   L6   0    0 ];
    dims = struct('L1', 0.333, 'L2', 0.316, 'L3', 0.0825, 'L4', 0.0825, 'L5', 0.384, 'L6', 0.088);
    mdh = mdh_table('Panda', rows, dims);
end
