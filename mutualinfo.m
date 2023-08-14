clc
clearvars
load pat-tau0.1.dat
C=reshape(pat_tau0_1,[128,301,128]);
n=128;
MI=[];
for snap=1:301
     D=C(:,snap,:);

    A=squeeze(D);
    for i=1:n
    for j=1:n
        if A(i,j)>2.5
           A(i,j)=1;
        else A(i,j)=0;
        end
    end
end
p1=sum(A,[1 2])/n^2;
p0=1-p1;
A1=[];
x=0;y=0;z=0;w=0; %%%%x1---p10 y1---p11 z1---p01 w1---p00
A1(1,1)=A(1,n);
A1(1,n+2)=A(1,1);
A1(1,2:n+1)=A(1,1:n);
A1(n+2,1)=A(1,n);
A1(n+2,n+2)=A(1,1);
A1(n+2,2:n+1)=A(1,1:n);

A1(1,1)=A(n,1);
A1(n+2,1)=A(1,1);
A1(2:n+1,1)=A(1:n,1);
A1(1,n+2)=A(n,1);
A1(n+2,n+2)=A(1,1);
A1(2:n+1,n+2)=A(1:n,1);
A1(2:n+1,2:n+1)=A;

for i=2:n+1
    for j=2:n+1
        if (A1(i,j)==1 && A1(i-1,j)==0);
                   x=x+1;
        end
              
        if (A1(i,j)==1 && A1(i+1,j)==0);
                   x=x+1;
          end
         if (A1(i,j)==1 && A1(i,j-1)==0);
                   x=x+1;
          end
           
          if (A1(i,j)==1 && A1(i,j+1)==0);
                   x=x+1;
          end
    end
end


for i=2:n+1
    for j=2:n+1
        if (A1(i,j)==1 && A1(i-1,j)==1);
                   y=y+1;
        end
              
        if (A1(i,j)==1 && A1(i+1,j)==1);
                   y=y+1;
          end
         if (A1(i,j)==1 && A1(i,j-1)==1);
                   y=y+1;
          end
           
          if (A1(i,j)==1 && A1(i,j+1)==1);
                   y=y+1;
          end
end
end
for i=2:n+1
    for j=2:n+1
        if (A1(i,j)==0 && A1(i-1,j)==1);
                   z=z+1;
        end
              
        if (A1(i,j)==0 && A1(i+1,j)==1);
                   z=z+1;
          end
         if (A1(i,j)==0 && A1(i,j-1)==1);
                   z=z+1;

          end
           
          if (A1(i,j)==0 && A1(i,j+1)==1);
                   z=z+1;

          end
end
end


for i=2:n+1
    for j=2:n+1
        if (A1(i,j)==0 && A1(i-1,j)==0);
                   w=w+1;
        end
              
        if (A1(i,j)==0 && A1(i+1,j)==0);
                    w=w+1;
          end
         if (A1(i,j)==0 && A1(i,j-1)==0);
                    w=w+1;

          end
           
          if (A1(i,j)==0 && A1(i,j+1)==0);
                   w=w+1;

          end
end
end

tot_pair=n^2*(n^2-1)/2;

x1=x/tot_pair;
y1=y/tot_pair;
z1=z/tot_pair;
w1=w/tot_pair;
%%%%%%%%%calculate mutul information%%%%%%%%%%%%
M=[];

for k=2:n+1
    for l=2:n+1
        pklb=[];pkij=[];plij=[];
        b1=A1(k-1,l);b2=A1(k+1,l);b3=A1(k,l-1);b4=A1(k,l+1);
        B=[b1,b2,b3,b4];  %%%neighbours
        if A1(k,l)==1
            
            for i1=1:4
                 pkij(i1)=p1;
                if B(i1)==1
                    pklb(i1)=y1;
                    plij(i1)=p1;
                else 
                    pklb(i1)=x1;
                    plij(i1)=p0;
                end
            end
        else
            
            for i1=1:4
                pkij(i1)=p0;
                if B(i1)==1
                    pklb(i1)=z1;
                    plij(i1)=p1;
                else 
                    pklb(i1)=w1;
                    plij(i1)=p0;
                end
            end
        end
        M(k,l)=sum(pklb(i1).*log(pklb./(pkij.*plij)));
    end
end
M1=sum(M,[1 2])./n^2;
MI=[MI;M1];
end
M1=(MI-min(MI))./(max(MI)-min(MI));   %%%%normalized MI
