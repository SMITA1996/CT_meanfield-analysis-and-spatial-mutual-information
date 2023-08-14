syms x y r a D N  xs u tau1 sigma1
q1=[];q11=[];q2=[];q3=[];
% a1=1:0.1:4;
x1=0.01:0.01:4;
for j=1:length(x1)
 a=1.5;
s=x1(end);
y=0.4286;     %%%%calculate y using spatial_bif_analytic.m for respective parameter values

sigma1=0.35;r=0.1;D=0.1;tau1=0.5;
% f1=@(x) r+((a.*x.^2)./(1+x.^2))-x ;
f=@(x) r+((a.*x.^2)./(1+x.^2))-x + (D.*(y-x));
Rexplicit = solve(f(x),x,'MaxDegree',3);
Rnumeric = vpa(Rexplicit);
xs=Rnumeric(1);

g(x)= diff(f(x));
u=1-tau1.*g(xs);

q4=@(x)(f(x)+(sigma1.*x./u))./(sigma1.*(x.^2)./u);
p1=inline(vectorize(q4(x)),'x');
q5=@(x)(f(x)./(sigma1.*(x.^2)./u));
q5=inline(vectorize(q5(x)),'x');
p2=@(x) (exp(int(p1(x),0.01,x)).*u)./(sigma1.*(x.^2));    %%%%%%%stationary pdf (p/q2)
p2=inline(vectorize(p2(x)),'x');
p=p2(x1(j));
%  p=double(p);
%  q2=integral(@(x) p2(x),0.01,s);
q2=integral(@(x) p2(x),0.01,4);  %%%normalization constant

q2=double(q2);

q5=@(x)(f(x)./(sigma1.*(x.^2)./u));
q6=@(x) int(q5,x,[0.04 x]);  %%% constant% q5=inline(vectorize(q5(x)),'x')

 phi=@(x)((0.5*log(sigma1*(x^2)./u))-q6(x));
 phi2=@(x) int(exp(-phi(x)),x,[0.04 x]);
phi=inline(vectorize(phi(x)),'x');
phi1=phi(s);
q1=[q1;   x1(j)  p./q2 phi1];
   j=length(x1);
end

