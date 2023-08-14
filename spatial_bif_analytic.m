clc
clear all
syms x y r a D N  xs u tau1 sigma1
q1=[];q11=[];q2=[];q3=[];
a1=3;        %% solving the below self consistent equation for a1=1:0.01:4 yields the analytical bifurcation diagram of the system
x1=0.01:0.01:4;
s=x1(end);
for ii=1:length(a1)
  ii
for j=1:length(x1)
y=x1(j);a=a1(ii);
sigma1=0.7;r=0.1;D=0.001;tau1=5;
% f1=@(x) r+((a.*x.^2)./(1+x.^2))-x ;
f=@(x) r+((a.*x.^2)./(1+x.^2))-x + (D.*(y-x));
Rexplicit = solve(f(x),x,'MaxDegree',3);
Rnumeric = vpa(Rexplicit);
xs=Rnumeric(1);
g(x)= diff(f(x));
u=1-tau1.*g(xs);

q4=@(x)(f(x)+(sigma1.*x./u))./(sigma1.*(x.^2)./u);
p1=inline(vectorize(q4(x)),'x');

 p2=@(x) (exp(int(p1(x),0.01,x)).*u)./(sigma1.*(x.^2));     %%%%%%%%%%% p2*Nc is Ps(x)
 p2=inline(vectorize(p2(x)),'x');


p3=@(x) x.*(exp(int(p1(x),0.01,x)).*u)./(sigma1.*(x.^2));
p3=inline(vectorize(p3(x)),'x');

  q2=integral(@(x) p2(x),0.01,s);
    q3=integral(@(x) p3(x),0.01,s);
    if abs(q3/q2-y)<0.01
         q11=[q11; x1(j)  q3/q2];

   j=length(x1);  
 end
   
end
end

