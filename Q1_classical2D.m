%% Monte Carlo simulation of a classical electron in a 2D oscillator
%%% Written by Ashwin Nayak, A53212870, asnayak[at]ucsd.edu
%% Initial Conditions
% Problem Parameters
dim = 2; beta = 1/200;

% Simulation Parameters
npar =10000; n_iter = 2000;
h = 5; % MC Step size

% Initial declarations
X = rand(dim,npar); P = rand(dim,npar);
E = zeros(n_iter,1); 
E(1) = avg_energy(X,P,npar);

%% Propagation
prog = waitbar(0,'Progress Bar');
tic
for i = 1:n_iter
    for ipar = 1:npar
        % Pick a random update
        dX=2*rand(dim,1)-1; dP = 2*rand(dim,1)-1;
    
        % Calculate acceptance probability
        rho_ratio = exp(-beta*(energy(X(:,ipar)+h*dX,P(:,ipar)+h*dP)-energy(X(:,ipar),P(:,ipar))));
        A = min(1,rho_ratio);
    
        % Update Position
        if (A>rand())
            X(:,ipar) = X(:,ipar)+h*dX; P(:,ipar) = P(:,ipar)+h*dP;
        end
    end
    E(i) = avg_energy(X,P,npar);
    prog = waitbar(i/n_iter,prog,['iteration: ',num2str(i)]);
end
toc
close(prog)     
%% Visualization
plot(E); hold on

%% Classical energy for a state
function E = energy(X,P)
    m = 1; w = 1;
    E = (0.5*((sum(P.^2)/m)+(m*(w^2)*sum(X.^2)))); 
    % E = E^2;
end

%% Average Energy of all configurations
function [E_avg] = avg_energy(X,P,npar)
    E_avg = 0; 
    for ipar=1:npar       
        E_avg = E_avg + energy(X(:,ipar),P(:,ipar));
    end
    E_avg = E_avg/npar; 
end