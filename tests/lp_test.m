function [pass] = lp_test()
    n = 2;
    m = 3;
    p = 0;
    P = [];
    c = [-1;-2];
    A = [];
    b = [];
    G = [-1 0;0 -1;1 1];
    h = [0;0;1];
    l = 3;
    nsoc = 0;
    q = [];
    solver = qoco;
    solver.setup(n, m, p, P, c, A, b, G, h, l, nsoc, q)
    out = solver.solve();

    tol = 1e-4;
    x_exp = [0.0;1.0];
    s_exp = [0.0;1.0;0.0];
    z_exp = [1.0;0.0;2.0];
    pass = (out.status == 1) && norm(out.x - x_exp) < tol && ...
        norm(out.s - s_exp) < tol && norm(out.z - z_exp) < tol;
end
