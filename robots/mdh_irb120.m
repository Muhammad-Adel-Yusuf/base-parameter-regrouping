function mdh = mdh_irb120()
%MDH_IRB120 Modified DH table of the ABB IRB 120 (q = 0 is the manufacturer's zero pose).
%   Columns: [alpha_{i-1} a_{i-1} d_i theta0_i].
    syms L1 L2 L3 L4 L5 real
    rows = [ 0      0    L1   0;
            -pi/2   0    0   -pi/2;
             0      L2   0    0;
            -pi/2   L3   L4   0;
             pi/2   0    0    0;
            -pi/2   0    L5   0 ];
    dims = struct('L1', 0.29, 'L2', 0.27, 'L3', 0.07, 'L4', 0.302, 'L5', 0.072);
    mdh = mdh_table('IRB120', rows, dims);
end
