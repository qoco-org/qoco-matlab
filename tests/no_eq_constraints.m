function [pass] = no_eq_constraints()
    n = 6;
    m = 6;
    p = 0;
    l = 3;
    nsoc = 1;
    P = diag([1,2,3,4,5,6]);
    c = [1;2;3;4;5;6];
    A = [];
    b = [];
    G = -eye(6);
    h = [0;0;0;0;0;0];
    q = [3];
    solver = qoco;
    solver.setup(n, m, p, P, c, A, b, G, h, l, nsoc, q)
    out = solver.solve();

    tol = 1e-3;
    x_exp = [0.0000;-0.0000;-0.0000;0.3981;-0.2624;-0.2993];
    s_exp = [0.0000;0.0000;0.0000;0.3981;-0.2624;-0.2993];
    z_exp = [1.0000;2.0000;3.0000;5.5923;3.6878;4.2040];
    pass = (out.status == 1) && norm(out.x - x_exp) < tol && ...
        norm(out.s - s_exp) < tol && ...
        norm(out.z - z_exp) < tol;
end
