function mdh = mdh_rpr_twist()
%MDH_RPR_TWIST Revolute-prismatic-revolute test robot with a GENERAL twist (neither parallel nor perpendicular).
%   Gautier and Khalil (1990), Section V-C, state that the parameters mX, mY, mZ of a translational link lying
%   between two revolute joints have a closed-form regrouping only if the joints in between are parallel or
%   perpendicular; otherwise supplementary regrouping has to be studied case by case. This table puts a
%   prismatic joint between two revolute ones with twists alpha = atan(3/4), whose sine and cosine (3/5, 4/5)
%   are rational, so that the exact certificates still apply. mdh_rpr_perp and mdh_rpr_par are the same robot
%   with perpendicular and parallel axes, for comparison. Columns: [alpha_{i-1} a_{i-1} d_i theta0_i type].
    syms L1 L2 L3 L4 L5 real
    tw = atan(sym(3)/4);
    rows = [ 0    0    L1   0   0;
             tw   L2   L3   0   1;
             tw   L4   L5   0   0 ];
    dims = struct('L1', 0.30, 'L2', 0.20, 'L3', 0.15, 'L4', 0.25, 'L5', 0.10);
    mdh = mdh_table('RPR_twist', rows, dims);
end
