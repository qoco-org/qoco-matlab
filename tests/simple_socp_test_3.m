function [pass] = simple_socp_test_3()
    n = 6;
    m = 6;
    p = 2;
    P = diag([1,2,3,4,5,6]);
    c = [1;2;3;4;5;6];
    A = [1 1 0 0 0 0;0 1 2 0 0 0];
    b = [1;2];
    G = -eye(6);
    h = [0;0;0;0;0;0];
    l = 0;
    nsoc = 2;
    q = [3,3];
    solver = qoco;
    solver.setup(n, m, p, P, c, A, b, G, h, l, nsoc, q)
    out = solver.solve();

    tol = 1e-3;
    x_exp = [1.0000;0.0000;1.0000;0.3981;-0.2625;-0.2993];
    s_exp = [1.0000;0.0000;1.0000;0.3981;-0.2625;-0.2993];
    y_exp = [4.0000;-6.0000];
    z_exp = [6.0000;-0.0000;-6.0000;5.5923;3.6876;4.2043];
    pass = (out.status == 1) && norm(out.x - x_exp) < tol && ...
        norm(out.s - s_exp) < tol && norm(out.y - y_exp) < tol && ...
        norm(out.z - z_exp) < tol;
end
