
% Paper to cite and code is based on
% [1]Kamchatnov, A., Kuo, Y., Lin, T., Horng, T., Gou, S., Clift, R., . . . Grimshaw, R. (2013).
% Transcritical flow of a stratified fluid over topography: 
% Analysis of the forced Gardner equation. 
% Journal of Fluid Mechanics, 736, 495-531. doi:10.1017/jfm.2013.556

% This code is Method of lines finite difference in space x and
% an ode solver is used to step forward in time

clear;
global nt h tau delta N up um alpha xend x0 epsilon gm l

tii=cputime;



%--------------------------------
%   Paramters for the problem
%
%------------------------------
% equation used -ut-delta*ux-6*alpha*u^2*ux+uxxx+Gx=0

tau=01.0;%time increments
alpha=-0.8% 4.8;%alpha
delta=1% -1;
gmx=(1)/(alpha^2)*(1-(2*alpha*delta)/(3))^(3/2) %eqn 2.12 [1]
gm=0.32%0.33;
l=10.0; % length of forcing l>>1
epsilon=1;
N=4000;% Number spatial pts
tend=85.0;




%spatial domain
xend=500;%xend
x0=-400;%x begin
h=abs((xend-x0)/N);%x increments
x=[x0:h:xend];% full x domain


%--------------------------------
%   Initial condition
%
%------------------------------
u=x;%quick way of making u have zeros say length as x
u=u*0;
u0=u;

%--------------------------------
%   Numerical solution
%------------------------------

plot(x,u0);%initial condition
tspan=[0:tau:tend];%time domain
for j=2:length(tspan)

    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    % Can play around with different matlab solvers
    [t un] = ode113('gkdvforce',[tspan(j-1) tspan(j)],u0,options);%finite difference scheme found in function gkdvforce
    u0=un(end,:);% put solution back into u0 to run again
   
    fprintf('time= %d \n ',t(end))
end

plot(x,un(end,:));
title(['Numerical Gardner    ',' t=  ',num2str(tend),'  \Delta= ',num2str(delta),' Gm= ',num2str(gm),' l= ',num2str(l),' alpha= ',num2str(alpha)])
timer = (cputime-tii)/60;
