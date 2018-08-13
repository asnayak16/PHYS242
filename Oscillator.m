%% Feynman Path Integral for the Harmonic Oscillator
%  Written by : Ashwin Nayak, asnayak@ucsd.edu
% -------------------------------------------------------------------------
%% Parameters
% Grid Parameters
Nd = 1000;
x0 = -4; xD = 4;
x = linspace(x0,xD,Nd+1);
dx = (xD - x0)/Nd;

% Time Parameters
T0 = 2*pi;
dt = T0/256; 
T_step = T0/256; 
N = T_step/dt;
T = T_step:T_step:5*T0;

% Wavefunction parameters
alp = 2; x_start = 0.75;
psi = NaN(Nd+1,17);
psi(:,1) = ((alp/pi)^0.25)*exp(-0.5*alp*(x-x_start).^2)';
psi(:,1) = psi(:,1)/(norm(psi(:,1))*sqrt(dx));
% -------------------------------------------------------------------------
%% Propagator Matrix 
A = sqrt(2*pi*1i*dt);
K_e = NaN(Nd+1,Nd+1);
for j=1:Nd+1
    K_e(j,:) = (1/A)*exp(1i*0.5*dt*(((x-x(j))/dt).^2 -((x+x(j))*0.5).^2));
end
K  = (dx^(N-1))*(K_e^N); 
% -------------------------------------------------------------------------
%% Wavefunction Propagation
% Preallocate Average Position Vector
avg_x = NaN(length(T)+1,1);
avg_x(1) = dx*x*(conj(psi(:,1)).*(psi(:,1)));

% Create Sparse Finite Difference Operator 
grad2 = spdiags([ones(Nd+1,1),-2*ones(Nd+1,1),ones(Nd+1,1)],[-1 0 1],Nd+1,Nd+1);
grad2(1,1) = 1; grad2(1,2) = -2; grad2(1,3) = 1;
grad2(end,end-2) = 1; grad2(end,end-1)=-2; grad2(end,end)=1;

% Preallocate Average energy arrays
KE = NaN(length(T)+1,1); PE = NaN(length(T)+1,1);
KE(1) = -(0.5/(dx))*(psi(:,1)')*grad2*psi(:,1);
PE(1) = 0.5*dx*(x.^2)*(conj(psi(:,1)).*psi(:,1));

fig1 = figure(1);
set(gcf,'Position',[100 -100 800 800]);
plt = plot(x,abs(psi(:,1)).^2,'r');%hold on;
set(gca,'XLim',[-4 4],'YLim',[0 0.9]); grid on; box on;
title('Time Evolution of Wavefunction in a Harmonic Oscillator');
xlabel('Position, x'); ylabel('Probability, |\psi|^2','Interpreter','Tex');
ax = gca; ax.FontName = 'Book Antiqua'; ax.FontSize = 12;
pos = [0.744 0.861225 0.147 0.0398];str = 'Time = 0.00 \pi'; 
txtbox = annotation('textbox',pos,'String',str,'FontSize',12,...
    'FontName','Book Antiqua','FitBoxToText','off');

vid = VideoWriter('Oscillator2.mp4','MPEG-4'); %vid.FrameRate=8; 
open(vid);

% Propagate wavefunction
for it = 1:length(T)
    
    psi(:,it+1) = dx * K * psi(:,it);
    psi(:,it+1) = psi(:,it+1)/(norm(psi(:,it+1))*sqrt(dx)); 
    
    % Calculate Average values 
    avg_x(it+1) = dx*x*(conj(psi(:,it+1)).*psi(:,it+1));
    KE(it+1) = -(psi(:,it+1)')*grad2*psi(:,it+1)/(2*dx);
    PE(it+1) = 0.5*dx*(x.^2)*(conj(psi(:,it+1)).*psi(:,it+1));
    
    % Update Plot
    %pause(0.5)
    plt.YData = abs(psi(:,it+1)).^2; 
    txtbox.String = ['Time = ', num2str(round(2*T(it)/T0,3),'%.2f'),' \pi'];
    drawnow ;
    
    % Get Frame and add to video
     ax = gca;
     ax.Units = 'pixels';
     pos = ax.Position;
     ti = ax.TightInset;
     rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
     image2 = getframe(gcf);
     writeVideo(vid,image2);
end
close(vid);
% -------------------------------------------------------------------------
%% Plots, Visualization and Video
% Plot Average Position
fig2 = figure(2);
plot([0 T],avg_x);%,[0 T],real(KE));
xlabel('Time, t'); ylabel('Average Position, <X>');
% Plot Average Energies
fig3 = figure(3);
plot([0 T],real(KE),[0 T],PE,[0,T],real(KE)+PE);
xlabel('Time, t'); ylabel('Energy');

