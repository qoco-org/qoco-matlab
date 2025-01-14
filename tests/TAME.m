function [pass] = TAME()
    n = 2;
    m = 2;
    p = 1;
    P = [2 -2;-2 2];
    c = [0;0];
    A = [1 1];
    b = [1];
    G = -eye(2);
    h = [0;0];
    l = 2;
    nsoc = 0;
    q = [];
    solver = qoco;
    solver.setup(n, m, p, P, c, A, b, G, h, l, nsoc, q)
    out = solver.solve();

    tol = 1e-4;
    x_exp = [0.5000;0.5000];
    obj_exp = 0;
    pass = (out.status == 1) && norm(out.x - x_exp) < tol && ...
        norm(out.obj - obj_exp) < tol;
end
