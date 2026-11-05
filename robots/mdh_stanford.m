function mdh = mdh_stanford()
%MDH_STANFORD Modified DH table of a Stanford-type arm (RRPRRR): two intersecting revolute shoulder
%   joints with an offset along axis 2, a prismatic arm joint, and a spherical wrist whose three axes
%   intersect. Columns: [alpha_{i-1} a_{i-1} d_i theta0_i type]; type 1 = prismatic (q3 is added to
%   d3 = L3). Written from the kinematic structure in the modified (proximal) convention; it is not
%   copied from a published table, so compare with one before quoting numbers from the literature.
    syms L1 L2 L3 L4 real
    rows = [ 0      0   L1   0   0;
             pi/2   0   L2   0   0;
            -pi/2   0   L3   0   1;
             0      0   L4   0   0;
             pi/2   0   0    0   0;
            -pi/2   0   0    0   0 ];
    dims = struct('L1', 0.412, 'L2', 0.154, 'L3', 0.20, 'L4', 0.10);
    mdh = mdh_table('Stanford', rows, dims);
end
