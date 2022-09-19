function udot=gkdvforce(time,u)
global h delta N alpha x0 xend epsilon gm l

epsilon=1;
udot=u;%udot is ut time derivative

N=length(u);
x=[x0:h:xend];


for i=4:N-3
         
  
  udot(i)=-delta*(1/(2*h))*(u(i+1)-u(i-1))+(1/(1*h))*(1/1)*(u(i+1)+u(i)+u(i-1))*(u(i+1)-u(i-1))+(1/(8*h^3))*(-u(i+3)+8*u(i+2)-13*u(i+1)+13*u(i-1)-8*u(i-2)+u(i-3))-alpha*(1/(1*h))*(1/1)*(1/3)*(u(i+1)+u(i)+u(i-1)).^2*(u(i+1)-u(i-1)) -2*x(i)*(1/l^2)*(gm)*exp(-x(i)^2/(l^2));
 
  %periodic conditions
  udot(1)=0;
  udot(2)=0;
            
  udot(end)=0;
  udot(end-1)=0;

end
     
