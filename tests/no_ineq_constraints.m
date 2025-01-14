function [pass] = no_ineq_constraints()
    n = 6;
    m = 0;
    p = 2;
    l = 0;
    nsoc = 0;
    P = diag([1,2,3,4,5,6]);
    c = [1;2;3;4;5;6];
    A = [1 1 0 0 0 0;0 1 2 0 0 0];
    b = [1;2];
    G = [];
    h = [];
    q = [];
    solver = qoco;
    solver.setup(n, m, p, P, c, A, b, G, h, l, nsoc, q)
    out = solver.solve();

    tol = 1e-4;
    x_exp = [0.2000;0.8000;0.6000;-1.0000;-1.0000;-1.0000];
    y_exp = [-1.200;-2.400];
    pass = (out.status == 1) && norm(out.x - x_exp) < tol && ...
        norm(out.y - y_exp) < tol;
end
