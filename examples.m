%% Generate Sobol inputs
rng(5);
n_samp=2^13;
n_dim=5;

p = sobolset(n_dim,'Skip',2^5);
p = scramble(p,'MatousekAffineOwen');
u= 2*net(p, n_samp)-1;
%% Initialise IPC object
ipc1=infoprocap.IPC(u,8);
%% Run photonic system with the input
