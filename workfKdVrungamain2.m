
% Paper to cite and code is based on
% [1]Kamchatnov, A., Kuo, Y., Lin, T., Horng, T., Gou, S., Clift, R., . . . Grimshaw, R. (2013).
% Transcritical flow of a stratified fluid over topography: 
% Analysis of the forced Gardner equation. 
% Journal of Fluid Mechanics, 736, 495-531. doi:10.1017/jfm.2013.556


% This code will anitmate the forced gardner equation just press play!!! :)



clear;
tii=cputime;
% N = 2^10+2^9;
N=1200;
x = linspace(-300,600,N); 
delta_x = x(2) - x(1);
delta_k = 2*pi/(N*delta_x);

k = [0:delta_k:N/2*delta_k,-(N/2-1)*delta_k:delta_k:-delta_k];%wavenumbers

 c=1;

% making initial conditionswith tanh's bit of aw rigmaroll
u = 1/2*c*(sech(sqrt(c)/2*(x+8))).^2;
up=1.0;
um=3.0;
a=up-um;
a2=um-up;
u1= 0.5*a*(1+tanh(x-400))+um;
u2= 0.5*a2*(1+tanh(x+400))+um;
u=(u1+u2);
u=(u1+u2)-min(u)+0.1;
u=u-2;
u=u*0;
plot(x,u)
alpha=0.8;
delta_t = 100000*0.4/N^2;
tmax =40.0; nmax = round(tmax/delta_t);
%-------------------------------------
del= -1.0;
U = fft(u);
gm=0.330;
l=2;
G=gm*exp(-x.^2/l^2);
G=fft(G);

%  Main bit This is a pseudo spectral method where the main equation is
%  transfered to fourier domain k which is then solved using a runga-kutta
%  method
%-ut-delta*ux-6*alpha*u^2*ux+uxxx+Gx=0
for n = 1:nmax
    t=n*delta_t
    % first we solve the linear part
    
    U = U.*exp(-1i*k.^3*delta_t);
   
    
    k1=delta_t*(3i*k.*fft(real(ifft(U)).^2))-delta_t*del*(1i*k.*fft(real(ifft(U))))+delta_t*(1i*k.*fft(real(ifft(G)))) ...
        -delta_t*(alpha*2i*k.*fft(real(ifft(U)).^3));
    k2=delta_t*(3i*k.*fft(real(ifft(U+0.5*k1)).^2))-delta_t*del*(1i*k.*fft(real(ifft(U+0.5*k1))))+delta_t*(1i*k.*fft(real(ifft(G)))) ...
        -delta_t*(alpha*2i*k.*fft(real(ifft(U+0.5*k1)).^3));
    k3=delta_t*(3i*k.*fft(real(ifft(U+0.5*k2)).^2))-delta_t*del*(1i*k.*fft(real(ifft(U+0.5*k2))))+delta_t*(1i*k.*fft(real(ifft(G)))) ...
        -delta_t*(alpha*2i*k.*fft(real(ifft(U+0.5*k2)).^3));
    k4=delta_t*(3i*k.*fft(real(ifft(U+k3)).^2))-delta_t*del*(1i*k.*fft(real(ifft(U+k3))))+delta_t*(1i*k.*fft(real(ifft(G)))) ...
        -delta_t*(alpha*2i*k.*fft(real(ifft(U+k3)).^3));
    
    U=U+(1/6)*(k1+2*k2+2*k3+k4);
  
    drawnow;
    plot(x,real(ifft(U)))
    title([' t= ',num2str(t),'L= ',num2str(l)])
    axis([-200 200 -1 1])
    %
    
end

plot(x,real(ifft(U)))
title([' t= ',num2str(t)])
% axis([-200 200 -0.6 0.6])

title(['Numerical KdV   ',' t=  ',num2str(tmax),'  \Delta= ',num2str(del),' Gm= ',num2str(gm),' l= ',num2str(l),' alpha= ',num2str(alpha)])
timer = (cputime-tii)/60

