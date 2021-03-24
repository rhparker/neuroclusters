f = 0.8;

% % no clusters
% a = 4;
% mee = 0.1;
% mei = mee;
% mie = -a*mee;
% mii = -a*mee;

% clusters 
% Nc = 4;
% mii = -1;
% mie = -1;
% mei = -((1-f)/f)*mii;
% mee = -Nc*((1-f)/f)*mii;

Nc = 4;
mee = 0.7*Nc;
mie = 0.7;
mii = -a*0.7;
mei = -a*0.7;

% specify N
N  = 20;
p  = floor( N*f/Nc );
ne = Nc*p;
ni = N-ne;
H = cluster(N, Nc, p, mee, mei, mie, mii); 
Jtilde = [ (p-1)*mee ni*mei ; ne*mie (ni-1)*mii ];

% % specify Nc, p, Nci
% Nc = 2;
% p = 10;
% Nci = 3;
% pi = round( (Nc*p/Nci)*(1-f)/f );
% ne = Nc*p;
% ni = Nci*pi;
% N = ne+ni;
% H = clusterEI(Nc, p, Nci, pi, mee, mei, mie, mii);    
    
l  = eig(H);
realvals = uniquetol(real(l));

gstar = sqrt(N)/realvals(end);
g = gstar + 0.5

l3 = (g/sqrt(N))*l-1;

%%
% timestepping

x0 = randn(N,1)*0.1;

fn = @(t, x)(-x + (H/sqrt(N))*tanh(g*x));

t = linspace(0,100,3000);
x = rk4(fn, x0, t);

J = Hx( H, g, x(:,end) );
lJ = eig(J);
realvalsJ = uniquetol(real(lJ));
xinh = abs( x(end,end) );
k = 1 - (gstar*xinh)^2;

figure;
subplot(2,2,1);
plot(real(l), imag(l), '.r', 'MarkerSize', 30);
subplot(2,2,2);
plot(real(l3), imag(l3), '.', 'MarkerSize', 30);
subplot(2,2,3);
plot(t, x(1:ne, :), '-b', t, x(ne+1:end,:), '-r');
subplot(2,2,4);
l4 = g*lJ/sqrt(N) - 1;
plot(real(l4),imag(l4), '.', 'MarkerSize', 30);



%% functions

% generate matrix H for linearization about origin
function H = cluster(N, Nc, p, mee, mei, mie, mii)
    if Nc*p >= N
        H = 0
    else
        ni = N - Nc*p;
        eblock = mee*clusterblock(Nc, p);        
        iblock = mii*clusterblock(1, ni);     
        eiblock = mei*ones( Nc*p, ni);
        ieblock = mie*ones( ni, Nc*p);  
        H = [ eblock eiblock ; ieblock iblock ];
    end
end

% allow excitatory and inhibitory clusters
function H = clusterEI(Nc, p, Nci, pi, mee, mei, mie, mii)
    ni = Nci*pi;
    ne = Nc*p;
    eblock = mee*clusterblock(Nc, p);        
    iblock = mii*clusterblock(Nci, pi);     
    eiblock = mei*ones( ne, ni);
    ieblock = mie*ones( ni, ne);  
    H = [ eblock eiblock ; ieblock iblock ];
end

% make matrix for n blocks of size p, will be np x np
function M = clusterblock(n, p)
    mask  = ones(p,p) - eye(p);
    block = ones(p,p) .* mask;
    bcell = repmat({block}, 1, n);
    M = blkdiag(bcell{:});
end

% matrix H(x) for arbitrary fixed point linearization
function M = Hx(H, g, x)
    row = sech(g*x').^2;
    M = H .* repmat( row, length(x), 1 );
end

% Runge-Kutta 4 ODE solver
% t is time grid
function u = rk4(f, u0, t)
    u = u0;
    h = t(2) - t(1);
    for index = 1:(length(t) - 1)
       k1 = h*f( t(index), u(:,end) );
       k2 = h*f( t(index)+h/2, u(:,end)+0.5*k1 );
       k3 = h*f( t(index)+h/2, u(:,end)+0.5*k2 ); 
       k4 = h*f( t(index)+h, u(:,end)+k3 );
       u = [ u  u(:,end)+(k1 + 2*k2 + 2*k3 + k4)/6 ];
    end
end