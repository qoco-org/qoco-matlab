function [pass] = no_constraints()
    n = 6;
    m = 0;
    p = 0;
    l = 0;
    nsoc = 0;
    P = diag([1,2,3,4,5,6]);
    c = [1;2;3;4;5;6];
    A = [];
    b = [];
    G = [];
    h = [];
    q = [];
    solver = qoco;
    solver.setup(n, m, p, P, c, A, b, G, h, l, nsoc, q)
    out = solver.solve();

    tol = 1e-3;
    x_exp = [-1, -1, -1, -1, -1, -1];
    pass = (out.status == 1) && norm(out.x - x_exp) < tol;
end
