%% Generate Sobol inputs
rng(5);

n_samp=2^12;
n_dim=2;
p = sobolset(n_dim,'Skip',2^10);
p = scramble(p,'MatousekAffineOwen');
u2= 2*net(p, n_samp)-1; % 2 dimensional sobol inputs

n_samp=2^13;
n_dim=5;
p = sobolset(n_dim,'Skip',2^10);
p = scramble(p,'MatousekAffineOwen');
u5= 2*net(p, n_samp)-1; % 5 dimensional sobol inputs

%% Initialise IPC object
ipc2=infoprocap.IPC(u2,14);   % ipc object for 2 dimensional inputs with max_deg=14
ipc5=infoprocap.IPC(u5,8);   % ipc object for 5 dimensional inputs with max_deg=8

%% Initialise photonic system
phot1=Phot_sys();
phot1.P_avg=6.6;  % average power of the laser
phot1.L=5;    % length of the fiber
phot1.updateParams();

%% Run photonic system
u_rep=repmat(u2,1,10);    % repeating features 10 times
phot1.run(u_rep);
X2=phot1.prepReadouts(5);  % measure readouts at 5m

u_rep=repmat(u5,1,4);    % repeating features 4 times
phot1.run(u_rep);
X5=phot1.prepReadouts(5);  % measure readouts at 5m

%% Calculate capacity with algorithm 1
[C2,~]=ipc2.estCap(X2,1);
[C5,~]=ipc5.estCap(X5,1);

%% Plot Capacity plots
Cmat=infoprocap.Plotter.capMat(ipc2,C2);   % Capacity matrix
Cd=infoprocap.Plotter.capBar(ipc5,C5);     % Capacity barplot





