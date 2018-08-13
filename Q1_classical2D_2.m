%% Monte Carlo simulation of a classical electron in a 2D oscillator
%%% Written by Ashwin Nayak, asnayak[at]ucsd.edu
%% Initial Conditions
% Problem Parameters
dim = 100; beta = 1;
w=0.15;
% Simulation Parameters
npar =1000; n_iter = 800;
h = 2; % MC Step size
% Initial declarations
X = 2*rand(dim,npar)-1; 
E = zeros(n_iter,1); 
E(1) = avg_energy(X,npar);
%% Propagation
prog = waitbar(0,'Progress Bar');
tic
succ=0; 
for i = 1:n_iter
    for ipar = 1:npar
        for d = 1:dim
            % Pick a random update
            temp = zeros(dim,1);
            temp(d) = h*(2*rand()-1);
             
            % Calculate acceptance probability
            rho_ratio = exp(-(energy(X(:,ipar)+temp)-energy(X(:,ipar))));
            A = min(1,rho_ratio);
             
            % Update Position
            if (A>rand())
                X(d,ipar) = X(d,ipar)+temp(d);
                succ= succ +1;
            end
        end
        E(i) = E(i)+mean( w^2 * (X(:,ipar).^2));
    end
        E(i) = E(i)/npar;
        prog = waitbar(i/n_iter,prog,['iteration: ',num2str(i)]);
end
Acceptance_Rate=(succ/ (n_iter*npar*dim))* 100
toc
close(prog)     
%% Visualization
plot(E); hold on
title('Ground State Energy');
xlabel('Iterations');
ylabel('Energy')
 
%% Energy for a state
function E = energy(X)
    kbT = 0.02; dim = 1000; w=0.15;
    dt = 1/(kbT*dim);
%     E = 0.5*(((X(1)-X(end))^2)/dt + dt*w^2*0.25*(X(1)+X(end))^2);
%     for i=2:length(X)
        E = E + 0.5*(((X(i)-X(i-1))^2)/dt + dt*w^2*0.25*((X(i)+X(i-1))^2));
%     end
    Y = circshift(X,1);
    E = sum(0.5*((((X-Y).^2)/dt) + (dt*(w^2)*0.25*((X+Y).^2))));
end
 
%% Average Energy of all configurations
function [E_avg] = avg_energy(X,npar)
    E_avg = 0; 
    for ipar=1:npar       
        E_avg = E_avg + energy(X(:,ipar));
    end
    E_avg = E_avg/npar; 
end
