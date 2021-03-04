% fraction excitatory
f = 0.8;
a = f/(1-f);

% clusters 
N  = 400;
Ne = round( f*N );
Ni = round( (1-f)*N );
Nc  = 1;
Nci = 10 ;

p  = round( Ne/Nc );
pi = round( Ni/Nci );

mm = 1;
mii = -mm;
mei = -mm;
mie = ((1-f)/f)*mm;
mee = ((1-f)/f)*mm;

H = clusterEI(Nc, p, Nci, pi, mee, mei, mie, mii);
Jtilde = [ (pi-1)*mii Ne*mie ; Ni*mei (Ne-1)*mee ];
Jeig = eig(Jtilde);
    
l  = eig(H);
realvals = uniquetol(real(l));

gstar = sqrt(N)/realvals(end-1);
g = gstar + 100;

l3 = (g/sqrt(N))*l-1;

% Jacobian manually
fn = @(t, x)(-x + (H/sqrt(N))*tanh(g*x));
dd = fn( 0, 0.001 * eye(N) ) / 0.001;
lcalc = eig(dd);

lambdac = (g/sqrt(N)) * (f*N/Nc - 1)*mee - 1;
lambdai = (g/sqrt(N)) * (-mii) - 1;
lambdar = (mee/2)*f*N*( (1/Nc - 1) + (a-1)/(f*N) );

%%
% timestepping

x0 = randn(N,1)*0.5;
% x0(1:16) = 0;
% x0(17:18) = 1;
% x0(19:20) = -1;

% x0(1:16) = 0.00088638 / 100;
% x0(17:19) = -0.0278 / 100;
% x0(20) = 0.09056 / 100;

fn = @(t, x)(-x + (H/sqrt(N))*tanh(g*x));

t = linspace(0,60,1000);
x = rk4(fn, x0, t);

J = Hx( H, g, x(:,end) );
lJ = eig(J);
realvalsJ = uniquetol(real(lJ));

xexc = abs( x(1,end) );
kexc = 1 - (gstar*xexc)^2;
xinh = abs( x(end,end) );
kinh = 1 - (gstar*xinh)^2;

lc = (g/sqrt(N))*(p-1)*kexc*mee - 1;
li = (g/sqrt(N))*a*mee - 1;

figure;
subplot(2,2,1);
plot(real(l), imag(l), '.r', 'MarkerSize', 30);
subplot(2,2,2);
plot(real(l3), imag(l3), '.', 'MarkerSize', 30);
subplot(2,2,3);
plot(t, x(1:f*N,:),'--b',t,x(f*N+1:end,:),'-r','LineWidth',1);
subplot(2,2,4);
l4 = g*lJ/sqrt(N) - 1;
plot(real(l4),imag(l4), '.', 'MarkerSize', 30);

%%

figure;
plot(t, x(1:f*N,:),'--b',t,x(f*N+1:end,:),'-r','LineWidth',1);
title('N = 200, Nci = 10, g = 10');
xlabel('t');


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