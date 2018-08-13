%% Monte Carlo simulation of a classical electron in a 2D oscillator
%%% Written by Ashwin Nayak, asnayak[at]ucsd.edu
%% Initial Conditions
clear; clc;

% Problem Parameters
dim = 100; kbT = 0.1; w=0.15;

% Simulation Parameters
npar =500; n_iter = 1000;
h = 1; % MC Step size

% Initial declarations
X = 2*rand(dim,npar)-1; 
E = zeros(n_iter,1);

%% Propagation
prog = waitbar(0,'Progress Bar');
tic
succ=0; total=0;
for i = 1:n_iter
    for ipar = 1:npar
        for d = 1:dim
            % Pick a random update
            temp = zeros(dim,1);
            temp(d) = h*(2*rand()-1);
            total=total+1;
             
            % Calculate acceptance probability
            rho_ratio = exp(-(energy(X(:,ipar)+temp,kbT,dim,w)-energy(X(:,ipar),kbT,dim,w)));
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
toc
close(prog) 

Acceptance_Rate=(succ/ total)* 100;

%% Calculating the Analytical probabilty Density
y=[-10:0.01:10];
l=sqrt(1/(w));
alpha=(w/pi)^0.25;
Analytical= alpha^2 * exp(-y.^2/l^2);

%% Calculate the Analytical Energy
kbT=100;
    for k=1:n_iter
        %only for the 
        E_analytical(k)=exp(-1/(w*kbT)) / (1-exp(-1/(w*kbT)));
    end

    
%% Visualization
figure(1);
plot(E);
title('Energy');
hold on;
plot(E_analytical);
legend('Simulation Energy','Analytical');
hold off;

figure(2);
histogram(X(:),100,'Normalization','pdf');
title('Histogram of X Distribution')
xlabel('Position');
ylabel('Frequency');
hold on;
plot(y,Analytical);
legend('Simulation Results','Analytical Results');
 
 
%% Energy for a state
function E = energy(X,kbT,dim,w)
    dt = 1/(kbT*dim);
    E = 0.5*(((X(1)-X(end))^2)/dt + dt*w^2*0.25*(X(1)+X(end))^2);
    for i=2:length(X)
        E = E + 0.5*(((X(i)-X(i-1))^2)/dt + dt*w^2*0.25*((X(i)+X(i-1))^2));
    end
%     Y = circshift(X,1);
%     E = sum(0.5*((((X-Y).^2)/dt) + (dt*(w^2)*0.25*((X+Y).^2))));
end
