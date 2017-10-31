%% Monte Carlo simulation : Anharmonic oscillator
%%% Written for PHYS 242 Final
%% clean up
clc; clear all;

%% Initial Conditions
% Problem parameters
a = 1; b = 1.4;       % Double Well parameters
kbT = 0.1;         % Temperature

% Simulation parameters
Npar  =  100        ;	% Number of Initial configurations  
Niter =  100       ;   % Number of MC iterations
Ndim  =  100        ;   % Number of points along path
h     =  1.2        ;	% MC Step size
dt = 1/(kbT*Ndim);

% Initial declarations
X = 2*rand(Ndim,Npar)-1; 
E = zeros(Niter,1); 
Naccept = 0;    % Acceptance Counter

%% Propagation
prog = waitbar(0,'Progress Bar');   % Progress bar
%profile on
tic
for iter = 1:Niter
    for ipar = 1:Npar
        for dim = 1:Ndim
            % Random update to one of the points on the path
            X_update = X(:,ipar);
            X_update(dim) = X_update(dim)+ h*(2*rand()-1);
             
            % Calculate acceptance probability
            rho_ratio = exp(-(energy(X_update,a,b,dt)-energy(X(:,ipar),a,b,dt)));
            A = min(1,rho_ratio);
             
            % Update Position
            if (A>rand())
                X(:,ipar) = X_update;
                Naccept = Naccept + 1;
            end
        end
        % Calculate energy for a path
        E(iter) = E(iter) + mean(a*((3*mean(X(:,ipar)).^4)-4*(b^2)*mean(X(:,ipar)).^2 + b^4));
    end
    E(iter)=E(iter)/Npar;   % Average out energy if multiple configs are used
    prog = waitbar(iter/Niter,prog,['iteration: ',num2str(iter)]);
end

%calculating the double well potential
y=[-7:0.1:7];
double_Well= a *( y.^2 - b^2).^2;

toc
close(prog)
%profile viewer
%% Results and Visualization
accept_rate = 100*Naccept/(Niter*Npar*Ndim);

figure(1);
plot(E);
title('Energy');

figure(2);
histogram(X(:),100,'Normalization','pdf');
title('Histogram of X Distribution')
xlabel('Position');
ylabel('Frequency');

figure(3);
plot(y,double_Well);
title('double well potential');
ylim([0 300]);

%% Function: Energy of a given Path
function E = energy(X,a,b,dt)
    Y = circshift(X,1);
    E = sum(0.5*(((X-Y).^2)./dt )+ ...
                  dt*a*(((X+Y)./2).^2.-(b^2)).^2);
end
