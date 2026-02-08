fprintf("**linear shooting method**\n");
%approximate a solution forboundary value problem  -y''=p(x)y'+q(x)y+r(x)
%with a<=x<=b, y(a)=alpha, y(b)=beta using N sub-intervals
a=input('enter a: ');
b=input('enter b: ');
alpha=input("enter alpha: ");
beta=input("enter beta: ");
N=input("enter N: ");


p=@(x) -2/x;
q=@(x) 2/x^2;
r=@(x) sin(log(x))/x^2;
%eq=diff(y,2)==p(x)*diff(y)+q(x)*y(x)+r(x);

h=(b-a)/N;

u=[alpha;0];
v=[0;1];

for i=0:N-1
    j=i+1; %j=1:N, u(a,0) does not exist in matlab
    x=a+(j-1)*h;

    k(1,1)=h*u(2,j);
    k(1,2)=h*(p(x)*u(2,j)+q(x)*u(1,j)+r(x));

    k(2,1)=h*(u(2,j)+1/2*k(1,2));
    k(2,2)=h*(p(x+h/2)*(u(2,j)+1/2*k(1,2))+q(x+h/2)*(u(1,j)+1/2*k(1,1))+r(x+h/2));

    k(3,1)=h*(u(2,j)+1/2*k(2,2));
    k(3,2)=h*(p(x+h/2)*(u(2,j)+1/2*k(2,2))+q(x+h/2)*(u(1,j)+1/2*k(2,1))+r(x+h/2));

    k(4,1)=h*(u(2,j)+k(3,2));
    k(4,2)=h*(p(x+h)*(u(2,j)+k(3,2))+q(x+h)*(u(1,j)+k(3,1))+r(x+h));

    u(1,j+1)=u(1,j)+1/6*(k(1,1)+2*k(2,1)+2*k(3,1)+k(4,1));
    u(2,j+1)=u(2,j)+1/6*(k(1,2)+2*k(2,2)+2*k(3,2)+k(4,2));
    

    dk(1,1)=h*v(2,j);
    dk(1,2)=h*(p(x)*v(2,j)+q(x)*v(1,j));

    dk(2,1)=h*(v(2,j)+1/2*dk(1,2));
    dk(2,2)=h*(p(x+h/2)*(v(2,j)+1/2*dk(1,2))+q(x+h/2)*(v(1,j)+1/2*dk(1,1)));
    
    dk(3,1)=h*(v(2,j)+1/2*dk(2,2));
    dk(3,2)=h*(p(x+h/2)*(v(2,j)+1/2*dk(2,2))+q(x+h/2)*(v(1,j)+1/2*dk(2,1)));
    
    dk(4,1)=h*(v(2,j)+dk(3,2));
    dk(4,2)=h*(p(x+h)*(v(2,j)+dk(3,2))+q(x+h)*(v(1,j)+dk(3,1)));
    
    v(1,j+1)=v(1,j)+1/6*(dk(1,1)+2*dk(2,1)+2*dk(3,1)+dk(4,1));
    v(2,j+1)=v(2,j)+1/6*(dk(1,2)+2*dk(2,2)+2*dk(3,2)+dk(4,2));

end

w(1,1)=alpha;
w(2,1)=(beta-u(1,N+1))/v(1,N+1);
disp([a,w(1,1),w(2,1)]);

for i=1:N
    j=i+1;
    W1=u(1,j)+w(2,1)*v(1,j);
    W2=u(2,j)+w(2,1)*v(2,j);
    x=a+(j-1)*h;

    w(1,j)=W1;
    w(2,j)=W2;
    disp([x,W1,W2]);
end