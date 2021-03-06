%% Monte Carlo simulation : Anharmonic oscillator
%%% Written for PHYS 242 Final

%% Initial Conditions
global a b kbT Npar Niter Ndim h dt
% Problem parameters
a = 2; b = 3;       % Double Well parameters
kbT = 0.01;         % Temperature

% Simulation parameters
Npar  =  100        ;	% Number of Initial configurations  
Niter =  100        ;   % Number of MC iterations
Ndim  =   74        ;   % Number of points along path
h     = 1.75        ;	% MC Step size
dt = 1/(kbT*Ndim);
% Initial declarations
X = 2*rand(Ndim,Npar)-1; 
E = zeros(Niter,1); 
Naccept = 0;    % Acceptance Counter

%% Propagation
prog = waitbar(0,'Progress Bar');   % Progress bar
profile on
tic
for iter = 1:Niter
    for ipar = 1:Npar
        for dim = 1:Ndim
            % Random update to one of the points on the path
            X_update = X(:,ipar);
            X_update(dim) = X_update(dim)+ h*(2*rand()-1);
             
            % Calculate acceptance probability
            rho_ratio = exp(-(energy(X_update)-energy(X(:,ipar))));
            A = min(1,rho_ratio);
             
            % Update Position
            if (A>rand())
                X(:,ipar) = X_update;
                Naccept = Naccept + 1;
            end
        end
        % Calculate energy for a path
        E(iter) = E(iter) + mean(a*((3*X(:,ipar).^4)-4*(b^2)*X(:,ipar).^2 + b^4));
    end
    E(iter)=E(iter)/Npar;   % Average out energy if multiple configs are used
    prog = waitbar(iter/Niter,prog,['iteration: ',num2str(iter)]);
end
toc
close(prog)
profile viewer
%% Results and Visualization
accept_rate = 100*Naccept/(Niter*Npar*Ndim);
plot(E)

%% Energy of a given Path
function E = energy(X)
    global a b dt
    Y = circshift(X,1);
    E = sum(0.5*(((X-Y).^2)./dt )+ ...
                  dt*a*(((X+Y)./2).^2.-(b^2)).^2);
end
