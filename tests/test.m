rng(10);

n = 20; % Number of variables
P = rand(n);
P = P*P';
c = rand(n, 1);

p = 5; % Number of equality constraints
A = rand(p, n);
b = rand(p, 1);

l = 5; % Dimension of LP cone
nsoc = 2; % Number of second order cones
q = [3, 4]; % Dimension of second order cones

m = l + sum(q);
G = rand(m, n);
h = rand(m, 1);

solver = qoco;
settings.verbose=1;
solver.setup(n, m, p, P, c, A, b, G, h, l, nsoc, q, settings)
solver.solve()
