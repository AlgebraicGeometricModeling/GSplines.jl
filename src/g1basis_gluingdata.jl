function loc_idx(m,k,i,j)
    return (k-1)*m^2+(j-1)*m+i
end

function umod(k,N)
    return mod(N+k-1,N)+1
end


function g1basis_EV_gluing(knt::Array,glu::Array,N::Int64)

m, deg = dim_deg(knt)

bas=spzeros(N*m^2,N+3);
#I fix b00=1 and all the other vertices to 0

a0=glu[1];
a1=glu[2];
a2=glu[3];

#B1
#b00

VAL=1;

for i in 1:N
    bas[loc_idx(m,i,1,1),1]=VAL; #b00=1
end

#b10

c=zeros(N); #Right hand of the system
c[1]=1; #Vector for the circulant matrix
c[N-1]=1;
c[N]=-a0;
T=Circulant(c);
T=convert(Matrix,T);
k=nullspace(T);
F=[T;k[:,1]';k[:,2]'];
RHS=zeros(N+2);
RHS[1:N].=(2-a0)*VAL;
#RHS[N+1]=-1.5;
#RHS[N+2]=-1.5;
sol=F\RHS;
push!(sol,sol[1]);
for i in 0:N-1
    bas[loc_idx(m,i+1,2,1),1]=sol[i+1]; #b10
    bas[loc_idx(m,i+1,1,2),1]=sol[i+2]; #b01
end

#b20

for i in 1:N
    bas[loc_idx(m,i,3,1),1]=1/2*bas[loc_idx(m,i,2,1),1]-VAL/10; #b20
    bas[loc_idx(m,i,1,3),1]=1/2*bas[loc_idx(m,i,1,2),1]-VAL/10; #b02
end

#b11

c=zeros(N); #Right hand of the system
c[1]=1; #Vector for the circulant matrix
c[2]=1;
T=Circulant(c);
T=convert(Matrix,T);
if N%2==0 #N even
    k=zeros(1,N); #It will be the kernel vector -1 1 -1 ... 1
    for i in 1:N
        k[i]=(-1)^i;
    end
    Z=[T;k];
    for i in 1:N
        c[i]=1/5*((a0-2*a1)*VAL)+(2-a0+2/5*a1)*bas[loc_idx(m,i,2,1),1]+4/5*a0*bas[loc_idx(m,i,3,1),1];
    end
    ff=[c;0];
    sol=Z\ff;

    for i in 1:N
        bas[loc_idx(m,i,2,2),1]=sol[i]; #b11
    end
else
    for i in 1:N
        c[i]=1/5*((a0-2*a1)*VAL)+(2-a0+2/5*a1)*bas[loc_idx(m,i,2,1),1]+4/5*a0*bas[loc_idx(m,i,3,1),1];
    end
    sol=T\c;

    for i in 1:N
        bas[loc_idx(m,i,2,2),1]=sol[i]; #b11
    end
end


#b21-b12 and b31-b13

for i in 1:N
    R=1/5*(-1/2*a0+a1-1/2*a2)*VAL+1/2*(a0-2*a1+1/5*a2)*bas[loc_idx(m,i,2,1),1]+(-a0+4/5*a1+2)*bas[loc_idx(m,i,3,1),1];
    bas[loc_idx(m,i,3,2),1]=1/2*R; #b21
    bas[loc_idx(m,umod(i-1,N),2,3),1]=1/2*R; #b12
    R1=-3/5*a2*bas[loc_idx(m,i,3,1),1];
    bas[loc_idx(m,i,4,2),1]=1/2*R1; #b31
    bas[loc_idx(m,umod(i-1,N),2,4),1]=1/2*R1; #b13
end



#B2-B3

#b00

for i in 1:N
    bas[loc_idx(m,i,1,1),2]=0;
    bas[loc_idx(m,i,1,1),3]=0;
end

#b10

v=zeros(N);
v[1]=1;
v[N-1]=1;
v[N]=-a0;

H=Circulant(v);
H=convert(Matrix,H);
sol=nullspace(H);
sol1=sol[:,1];
sol1=sol1./maximum(sol1);
#push!(sol1,sol1[1]);
sol2=sol[:,2];
sol2=sol2./maximum(sol2);
#push!(sol2,sol2[1]);

for i in 1:N
    bas[loc_idx(m,i,2,1),2]=sol1[i];
    bas[loc_idx(m,umod(i-1,N),1,2),2]=sol1[i];
    bas[loc_idx(m,i,2,1),3]=sol2[i];
    bas[loc_idx(m,umod(i-1,N),1,2),3]=sol2[i];
    #=bas[i*36+2,2]=sol1[i+1]; #b10
    bas[i*36+7,2]=sol1[i+2]; #b01
    bas[i*36+2,3]=sol2[i+1]; #b10
    bas[i*36+7,3]=sol2[i+2]; #b01=#
end

#b20

for i in 1:N
    bas[loc_idx(m,i,3,1),2]=1/2*bas[loc_idx(m,i,2,1),2];
    bas[loc_idx(m,i,1,3),2]=1/2*bas[loc_idx(m,i,1,2),2];
    bas[loc_idx(m,i,3,1),3]=1/2*bas[loc_idx(m,i,2,1),3];
    bas[loc_idx(m,i,1,3),3]=1/2*bas[loc_idx(m,i,1,2),3];
    #=bas[i*36+3,2]=1/2*bas[i*36+2,2]; #b20
    bas[i*36+13,2]=1/2*bas[i*36+7,2]; #b02
    bas[i*36+3,3]=1/2*bas[i*36+2,3]; #b20
    bas[i*36+13,3]=1/2*bas[i*36+7,3]; #b02=#
end



#b11

c=zeros(N); #Right hand of the system
c1=zeros(N);
c2=zeros(N);
c[1]=1; #Vector for the circulant matrix
c[2]=1;
T=Circulant(c);
T=convert(Matrix,T);
if N%2==0 #N even
    k=zeros(1,N); #It will be the kernel vector -1 1 -1 ... 1
    for i in 1:N
        k[i]=(-1)^i;
    end
    Z=[T;k];
    for i in 1:N
        c1[i]=(2-a0+2/5*a1)*bas[loc_idx(m,i,2,1),2]+4/5*a0*bas[loc_idx(m,i,3,1),2];
        c2[i]=(2-a0+2/5*a1)*bas[loc_idx(m,i,2,1),3]+4/5*a0*bas[loc_idx(m,i,3,1),3];
        #=c1[i+1]=(1/5)*(-5*(a-2)*bas[i*36+2,2]+4*a*bas[i*36+3,2]);
        c2[i+1]=(1/5)*(-5*(a-2)*bas[i*36+2,3]+4*a*bas[i*36+3,3]);=#
    end
    cc1=[c1;0];
    cc2=[c2;0];
    sol1_11=Z\cc1;
    sol2_11=Z\cc2;

    for i in 1:N
        bas[loc_idx(m,i,2,2),2]=sol1_11[i];
        bas[loc_idx(m,i,2,2),3]=sol2_11[i]; #b11
    end

else

    for i in 1:N
        c1[i]=(2-a0+2/5*a1)*bas[loc_idx(m,i,2,1),2]+4/5*a0*bas[loc_idx(m,i,3,1),2];
        c2[i]=(2-a0+2/5*a1)*bas[loc_idx(m,i,2,1),3]+4/5*a0*bas[loc_idx(m,i,3,1),3];
    end
    sol1_11=T\c1;
    sol2_11=T\c2;

    for i in 1:N
        bas[loc_idx(m,i,2,2),2]=sol1_11[i];
        bas[loc_idx(m,i,2,2),3]=sol2_11[i]; #b11
    end
end


#b21-b12 and b31-b13


for i in 1:N

    R=1/2*(a0-2*a1+1/5*a2)*bas[loc_idx(m,i,2,1),2]+(-a0+4/5*a1+2)*bas[loc_idx(m,i,3,1),2];
    bas[loc_idx(m,i,3,2),2]=1/2*R; #b21
    bas[loc_idx(m,umod(i-1,N),2,3),2]=1/2*R; #b12
    R1=-3/5*a2*bas[loc_idx(m,i,3,1),2];
    bas[loc_idx(m,i,4,2),2]=1/2*R1; #b31
    bas[loc_idx(m,umod(i-1,N),2,4),2]=1/2*R1; #b13

    R=1/2*(a0-2*a1+1/5*a2)*bas[loc_idx(m,i,2,1),3]+(-a0+4/5*a1+2)*bas[loc_idx(m,i,3,1),3];
    bas[loc_idx(m,i,3,2),3]=1/2*R; #b21
    bas[loc_idx(m,umod(i-1,N),2,3),3]=1/2*R; #b12
    R1=-3/5*a2*bas[loc_idx(m,i,3,1),3];
    bas[loc_idx(m,i,4,2),3]=1/2*R1; #b31
    bas[loc_idx(m,umod(i-1,N),2,4),3]=1/2*R1; #b13

end



#B4-B6

#b11

for i in 0:N-1

    bas[loc_idx(m,i+1,2,2),4+i]=VAL; #b11^0=1

    bas[loc_idx(m,i+1,3,1),4+i]=(5/(4*a0))*bas[loc_idx(m,i+1,2,2),4+i]*VAL; #b20
    bas[loc_idx(m,umod(i,N),1,3),4+i]=(5/(4*a0))*bas[loc_idx(m,i+1,2,2),4+i]*VAL; #b02

    bas[loc_idx(m,i+1,1,3),4+i]=(5/(4*a0))*bas[loc_idx(m,i+1,2,2),4+i]*VAL; #b20
    bas[loc_idx(m,(i+1)%N+1,3,1),4+i]=(5/(4*a0))*bas[loc_idx(m,i+1,2,2),4+i]*VAL; #b02


#b30

bas[loc_idx(m,i+1,4,1),4+i]=bas[loc_idx(m,i+1,3,1),4+i]; #b30
bas[loc_idx(m,umod(i,N),1,4),4+i]=bas[loc_idx(m,i+1,3,1),4+i]; #b03

bas[loc_idx(m,i+1,1,4),4+i]=bas[loc_idx(m,i+1,3,1),4+i]; #b30
bas[loc_idx(m,(i+1)%N+1,4,1),4+i]=bas[loc_idx(m,i+1,3,1),4+i]; #b03


#b21-b12 b31-b13

for k in 0:1

    pidx=i+1+k;
    if pidx==N+1
        pidx=1;
    end

    R=(-a0+4/5*a1+2)*bas[loc_idx(m,pidx,3,1),4+i]+3/5*a0*bas[loc_idx(m,pidx,4,1),4+i];
    #R=-(a-2)*bas[(i+k)%N*36+3,4+i]+(3/5)*a*bas[i*36+4,4+i];
    bas[loc_idx(m,pidx,3,2),4+i]=1/2*R; #b21
    bas[loc_idx(m,umod(pidx-1,N),2,3),4+i]=1/2*R; #b12

    #R=3/5*a0*bas[loc_idx(m,pidx,3,1),4+i];
    #R=-(a-2)*bas[(i+k)%N*36+3,4+i]+(3/5)*a*bas[i*36+4,4+i];
    R=(a0-2*a1+2/5*a2)*bas[loc_idx(m,pidx,3,1),4+i]+(-a0+6/5*a1+2)*bas[loc_idx(m,pidx,4,1),4+i];
    bas[loc_idx(m,pidx,4,2),4+i]=1/2*R; #b21
    bas[loc_idx(m,umod(pidx-1,N),2,4),4+i]=1/2*R; #b12

    #=R=a*bas[(i+k)%N*36+3,4+i]-(a-2)*bas[i*36+4,4+i];
    bas[(i+k)%N*36+10,4+i]=1/2*R; #b31
    bas[(N-1+k+i)%N*36+20,4+i]=1/2*R; #b13=#
end

end

#=
for i in 0:N-1
        bas[i%N*36+10,7+i]=1; #b21
        bas[i%N*36+4,7+i]=1/2; #b21
        bas[(N-1+i)%N*36+19,7+i]=1/2; #b12
end
=#

return bas;

end

function g1basis_borderRV_gluing(knt::Array,glu::Array)

m, deg = dim_deg(knt)

bas=spzeros(2*m^2,4)

bas[loc_idx(m,1,1,1),1]=1/2
bas[loc_idx(m,2,1,1),1]=1/2
bas[loc_idx(m,1,2,1),1]=1

bas[loc_idx(m,1,1,2),2]=1/2
bas[loc_idx(m,2,2,1),2]=1/2
bas[loc_idx(m,1,2,2),2]=1

bas[loc_idx(m,1,1,2),3]=1/2
bas[loc_idx(m,2,2,1),3]=1/2
bas[loc_idx(m,2,2,2),3]=1

bas[loc_idx(m,1,1,1),4]=1/2
bas[loc_idx(m,2,1,1),4]=1/2
bas[loc_idx(m,2,1,2),4]=1


a0=glu[1];
a1=glu[2];

for k in 1:4
    bas[loc_idx(m,1,1,3),k]=1/2*bas[loc_idx(m,1,1,2),k]-1/10*bas[loc_idx(m,1,1,1),k]; #b20
    bas[loc_idx(m,2,3,1),k]=bas[loc_idx(m,1,1,3),k]; #b02
    bas[loc_idx(m,1,2,3),k]=1/2*((-a0+6/5*a1+2)*bas[loc_idx(m,1,1,3),k]+2/5*a0*bas[loc_idx(m,1,1,2),k]); #b21
    bas[loc_idx(m,2,3,2),k]=bas[loc_idx(m,1,2,3),k]; #b12
    bas[loc_idx(m,1,2,4),k]=3/10*a0*bas[loc_idx(m,1,1,3),k]; #b31
    bas[loc_idx(m,2,4,2),k]=bas[loc_idx(m,1,2,4),k]; #b13
end

return bas

end


function g1basis_noncrossingRV(kn::Array, gluL::Array, gluR::Array, i::Int64, j::Int64, N::Int64)

    m, d = dim_deg(kn);
    nbas=4;
    bas=spzeros(N*m^2,nbas*N);



    a=gluR[1];
    a2R=gluR[2];
    a3R=gluR[3];
    a2L=-gluL[2];
    a3L=-gluL[1];

    #=a=gluR[1];
    a2R=gluR[2];
    a3R=gluR[3];
    a2L=gluL[2];
    a3L=gluL[1];=#

    #B1
    #b00=1
    P=0; # value for b30(1)=b30(3)
    for k in 1:N


        bas[loc_idx(m,k,i+3,j),1+nbas*(k-1)]=P; #b31_1
        bas[loc_idx(m,umod(k-1,N),j,i+3),1+nbas*(k-1)]=P;

        bas[loc_idx(m,k,i-3,j),1+nbas*(k-1)]=P; #b31_3
        bas[loc_idx(m,umod(k-1,N),j,i-3),1+nbas*(k-1)]=P;



        bas[loc_idx(m,k,i,j),1+nbas*(k-1)]=1;
        bas[loc_idx(m,umod(k-1,N),j,i),1+nbas*(k-1)]=1;

        #b10
        bas[loc_idx(m,k,i+1,j),1+nbas*(k-1)]=1;
        bas[loc_idx(m,umod(k-1,N),j,i+1),1+nbas*(k-1)]=1;

        bas[loc_idx(m,k,i,j+1),1+nbas*(k-1)]=1;

        bas[loc_idx(m,k,i-1,j),1+nbas*(k-1)]=1;
        bas[loc_idx(m,umod(k-1,N),j,i-1),1+nbas*(k-1)]=1;

        bas[loc_idx(m,umod(k-1,N),j+1,i),1+nbas*(k-1)]=1;

        #b20
        bas[loc_idx(m,k,i+2,j),1+nbas*(k-1)]=2/5+P;
        bas[loc_idx(m,umod(k-1,N),j,i+2),1+nbas*(k-1)]=2/5+P;
        bas[loc_idx(m,k,i-2,j),1+nbas*(k-1)]=2/5+P;
        bas[loc_idx(m,umod(k-1,N),j,i-2),1+nbas*(k-1)]=2/5+P;
        #bas[loc_idx(m,k,i,j+2),1+3*(k-1)]=2/5;
        #bas[loc_idx(m,umod(k-1,N),j+2,i),1+3*(k-1)]=2/5;

        #b11
        bas[loc_idx(m,k,i+1,j+1),1+nbas*(k-1)]=1+a/25*(10*P-6); #b11_1
        bas[loc_idx(m,k,i-1,j+1),1+nbas*(k-1)]=1+a/25*(-10*P+6); #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),1+nbas*(k-1)]=1+a/25*(-10*P+6); #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),1+nbas*(k-1)]=1+a/25*(10*P-6); #b11_0


        #bas[loc_idx(m,k,i+1,j+1),1+nbas*(k-1)]=1-6/25*a; #b11_1
        #bas[loc_idx(m,k,i-1,j+1),1+nbas*(k-1)]=1+6/25*a; #b11_2
        #bas[loc_idx(m,umod(k-1,N),j+1,i-1),1+nbas*(k-1)]=1+6/25*a; #b11_3
        #bas[loc_idx(m,umod(k-1,N),j+1,i+1),1+nbas*(k-1)]=1-6/25*a; #b11_0




        #b21_12
        #bR=2/5*(1-3/5*a2R);
        #bL=2/5*(1-3/5*a2L);
        #bzero=2/5;

        bR=(1/2)*(1/25*((-10*a+20*a2R+50)*P-12*a2R)+4/5);
        bL=(1/2)*(1/25*((20*a2L+10*a+50)*P-12*a2L)+4/5);

        bas[loc_idx(m,k,i+2,j+1),1+nbas*(k-1)]=bR; #b21_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+2),1+nbas*(k-1)]=bR;

        #bas[loc_idx(m,k,i+1,j+2),1+3*(k-1)]=bzero; #b21_2
        #bas[loc_idx(m,k,i-1,j+2),1+3*(k-1)]=bzero;

        bas[loc_idx(m,k,i-2,j+1),1+nbas*(k-1)]=bL; #b21_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-2),1+nbas*(k-1)]=bL;

        #bas[loc_idx(m,umod(k-1,N),j+2,i+1),1+3*(k-1)]=bzero; #b21_0
        #bas[loc_idx(m,umod(k-1,N),j+2,i-1),1+3*(k-1)]=bzero;

        #b31_13
        #bbR=-6/25*a3R;
        #bbL=-6/25*a3L;

        bbR=(1/2)*(1/25*((-20*a2R+10*a3R+50)*P+6*a3R));
        bbL=(1/2)*(1/25*((-20*a2L+10*a3L+50)*P-6*a3L));

        bas[loc_idx(m,k,i+3,j+1),1+nbas*(k-1)]=bbR; #b31_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+3),1+nbas*(k-1)]=bbR;

        bas[loc_idx(m,k,i-3,j+1),1+nbas*(k-1)]=bbL; #b31_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-3),1+nbas*(k-1)]=bbL;

        #B2
        #b00=0

        bas[loc_idx(m,k,i,j),2+nbas*(k-1)]=0;
        bas[loc_idx(m,umod(k-1,N),j,i),2+nbas*(k-1)]=0;

        #b10
        #bas[loc_idx(m,k,i+1,j),2+3*(k-1)]=0; #b10_1
        #bas[loc_idx(m,umod(k-1,N),j,i+1),2+3*(k-1)]=0;

        bas[loc_idx(m,k,i,j+1),2+nbas*(k-1)]=-1; #b10_2

        #bas[loc_idx(m,k,i-1,j),2+3*(k-1)]=0; #b10_3
        #bas[loc_idx(m,umod(k-1,N),j,i-1),2+3*(k-1)]=0;

        bas[loc_idx(m,umod(k-1,N),j+1,i),2+nbas*(k-1)]=1; #b10_0

        #b20
        #bas[loc_idx(m,k,i+2,j),2+3*(k-1)]=-1/2; #b20_1
        #bas[loc_idx(m,umod(k-1,N),j,i+2),2+3*(k-1)]=-1/2;
        #bas[loc_idx(m,k,i-2,j),2+3*(k-1)]=1/2; #b20_3
        #bas[loc_idx(m,umod(k-1,N),j,i-2),2+3*(k-1)]=1/2;
        #bas[loc_idx(m,k,i,j+2),2+3*(k-1)]=0; #b20_2
        #bas[loc_idx(m,umod(k-1,N),j+2,i),2+3*(k-1)]=0; #b20_0

        #b11
        bas[loc_idx(m,k,i+1,j+1),2+nbas*(k-1)]=-1; #b11_1
        bas[loc_idx(m,k,i-1,j+1),2+nbas*(k-1)]=-1; #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),2+nbas*(k-1)]=1; #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),2+nbas*(k-1)]=1; #b11_0

        #b21_12
        #bR=1/2*(-1+1/5*(3*a2R-1/2*a3R));
        #bL=1/2*(1-1/5*(3*a2L-1/2*a3L));

        bas[loc_idx(m,k,i+2,j+1),2+nbas*(k-1)]=0; #b21_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+2),2+nbas*(k-1)]=0; #-1

        bas[loc_idx(m,k,i-2,j+1),2+nbas*(k-1)]=0; #b21_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-2),2+nbas*(k-1)]=0;

        #b31_13
        #bbR=3/10*a3R;
        #bbL=-3/10*a3L;

        bas[loc_idx(m,k,i+3,j+1),2+nbas*(k-1)]=0; #b31_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+3),2+nbas*(k-1)]=0;

        bas[loc_idx(m,k,i-3,j+1),2+nbas*(k-1)]=0; #b31_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-3),2+nbas*(k-1)]=0;

        #B3
        #b11
        bas[loc_idx(m,k,i+1,j+1),3+nbas*(k-1)]=1; #b11_1
        bas[loc_idx(m,k,i-1,j+1),3+nbas*(k-1)]=-1; #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),3+nbas*(k-1)]=1; #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),3+nbas*(k-1)]=-1; #b11_0

        #bas[loc_idx(m,k,i+2,j+1),3+nbas*(k-1)]=1; #b21_1
        #bas[loc_idx(m,umod(k-1,N),j+1,i+2),3+nbas*(k-1)]=-1;

        #bas[loc_idx(m,k,i-2,j+1),3+nbas*(k-1)]=-1; #b21_3
        #bas[loc_idx(m,umod(k-1,N),j+1,i-2),3+nbas*(k-1)]=1;


        #bas[loc_idx(m,k,i+3,j+1),3+nbas*(k-1)]=1; #b31_1
        #bas[loc_idx(m,umod(k-1,N),j+1,i+3),3+nbas*(k-1)]=-1;

        #bas[loc_idx(m,k,i-3,j+1),3+nbas*(k-1)]=-1; #b31_3
        #bas[loc_idx(m,umod(k-1,N),j+1,i-3),3+nbas*(k-1)]=1;


        #B4 Here b30=(T 0 0 0) , T parameter

        T=1;

        bas[loc_idx(m,k,i+3,j),4+nbas*(k-1)]=T;
        bas[loc_idx(m,umod(k-1,N),j,i+3),4+nbas*(k-1)]=T;

        #bas[loc_idx(m,k,i-3,j),4+nbas*(k-1)]=-1;
        #bas[loc_idx(m,umod(k-1,N),j,i-3),4+nbas*(k-1)]=-1;


        f=(a*T)/(a2L+8*a-a2R);

        #b10
        bas[loc_idx(m,k,i+1,j),4+nbas*(k-1)]=f*2;
        bas[loc_idx(m,umod(k-1,N),j,i+1),4+nbas*(k-1)]=f*2;

        bas[loc_idx(m,k,i,j+1),4+nbas*(k-1)]=f*a;

        bas[loc_idx(m,k,i-1,j),4+nbas*(k-1)]=-2*f;
        bas[loc_idx(m,umod(k-1,N),j,i-1),4+nbas*(k-1)]=-2*f;

        bas[loc_idx(m,umod(k-1,N),j+1,i),4+nbas*(k-1)]=f*a;

        #b20
        bas[loc_idx(m,k,i+2,j),4+nbas*(k-1)]=f+T;
        bas[loc_idx(m,umod(k-1,N),j,i+2),4+nbas*(k-1)]=f+T;
        bas[loc_idx(m,k,i-2,j),4+nbas*(k-1)]=-f;
        bas[loc_idx(m,umod(k-1,N),j,i-2),4+nbas*(k-1)]=-f;

        g=(a*T)/(40*a+5*(a2L-a2R));
        #b11
        bas[loc_idx(m,k,i+1,j+1),4+nbas*(k-1)]=g*(13*a+10+2*a2L); #b11_1
        bas[loc_idx(m,k,i-1,j+1),4+nbas*(k-1)]=-g*(3*a+2*a2L+10); #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),4+nbas*(k-1)]=-g*(3*a+2*a2L+10); #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),4+nbas*(k-1)]=g*(13*a+10+2*a2L); #b11_0

        #b21_12

        h=T/(40*a+5*(a2L-a2R));

        bR=(h/2)*(-2*(8*a^2+(a2L-14*a2R-(a3R/2)-45)*a-2*(a2R+(5/2))*(a2L-a2R)));
        bL=(h/2)*(a*(6*a2L-a3L-10));

        bas[loc_idx(m,k,i+2,j+1),4+nbas*(k-1)]=bR; #b21_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+2),4+nbas*(k-1)]=bR;

        bas[loc_idx(m,k,i-2,j+1),4+nbas*(k-1)]=bL; #b21_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-2),4+nbas*(k-1)]=bL;

        #b31_13

        bbR=(h/2)*(-4*(-a2R^2+(a2L+8*a2R+(a3R/2)+(5/2))*a2R+((-a3R/2)-(5/2))*a2L+a*(-13*(a3R/4)-20)));
        bbL=(h/2)*(3*a*a3L);

        bas[loc_idx(m,k,i+3,j+1),4+nbas*(k-1)]=bbR; #b31_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+3),4+nbas*(k-1)]=bbR;

        bas[loc_idx(m,k,i-3,j+1),4+nbas*(k-1)]=bbL; #b31_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-3),4+nbas*(k-1)]=bbL;

        #B4  Here b30=( 1 0 -1 0 )

        #=bas[loc_idx(m,k,i+3,j),4+nbas*(k-1)]=1;
        bas[loc_idx(m,umod(k-1,N),j,i+3),4+nbas*(k-1)]=1;

        bas[loc_idx(m,k,i-3,j),4+nbas*(k-1)]=-1;
        bas[loc_idx(m,umod(k-1,N),j,i-3),4+nbas*(k-1)]=-1;


        f=2*a/(8*a+a2L-a2R);

        #b10
        bas[loc_idx(m,k,i+1,j),4+nbas*(k-1)]=f*2;
        bas[loc_idx(m,umod(k-1,N),j,i+1),4+nbas*(k-1)]=f*2;

        bas[loc_idx(m,k,i,j+1),4+nbas*(k-1)]=f*a;

        bas[loc_idx(m,k,i-1,j),4+nbas*(k-1)]=-2*f;
        bas[loc_idx(m,umod(k-1,N),j,i-1),4+nbas*(k-1)]=-2*f;

        bas[loc_idx(m,umod(k-1,N),j+1,i),4+nbas*(k-1)]=f*a;

        #b20
        bas[loc_idx(m,k,i+2,j),4+nbas*(k-1)]=f+1;
        bas[loc_idx(m,umod(k-1,N),j,i+2),4+nbas*(k-1)]=f+1;
        bas[loc_idx(m,k,i-2,j),4+nbas*(k-1)]=-(f+1);
        bas[loc_idx(m,umod(k-1,N),j,i-2),4+nbas*(k-1)]=-(f+1);

        g=2*a/(40*a+5*(a2L-a2R));
        #b11
        bas[loc_idx(m,k,i+1,j+1),4+nbas*(k-1)]=g*(5*a+a2R+a2L+10); #b11_1
        bas[loc_idx(m,k,i-1,j+1),4+nbas*(k-1)]=-g*(-5*a+a2R+a2L+10); #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),4+nbas*(k-1)]=-g*(-5*a+a2R+a2L+10); #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),4+nbas*(k-1)]=g*(5*a+a2R+a2L+10); #b11_0

        #b21_12

        h=40*a+5*(a2L-a2R);

        bR=1/(2*h)*(-16*a^2+(-2*a2L+22*a2R+2*a3R+100)*a+4*(a2L-a2R)*(a2R+5/2));
        bL=1/(2*h)*(-16*a^2+(-22*a2L-2*a3L+2*a2R-100)*a-4*(a2L-a2R)*(a2L+5/2));

        bas[loc_idx(m,k,i+2,j+1),4+nbas*(k-1)]=bR; #b21_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+2),4+nbas*(k-1)]=bR;

        bas[loc_idx(m,k,i-2,j+1),4+nbas*(k-1)]=bL; #b21_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-2),4+nbas*(k-1)]=bL;

        #b31_13

        bbR=1/(2*h)*(4*a2R^2+(-4*a2L-32*a-2*a3R-10)*a2R+(2*a3R+10)*a2L+10*a*(a3R+8));
        bbL=1/(2*h)*(4*a2L^2+(-2*a3L+32*a-4*a2R-10)*a2L+(-10*a3L-80)*a+2*a2R*(a3L+5));

        bas[loc_idx(m,k,i+3,j+1),4+nbas*(k-1)]=bbR; #b31_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+3),4+nbas*(k-1)]=bbR;

        bas[loc_idx(m,k,i-3,j+1),4+nbas*(k-1)]=bbL; #b31_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-3),4+nbas*(k-1)]=bbL;=#


        #=#B5  Here b30=( 0 0 1 0 )

        bas[loc_idx(m,k,i-3,j),5+nbas*(k-1)]=1;
        bas[loc_idx(m,umod(k-1,N),j,i-3),5+nbas*(k-1)]=1;

        f=a/(8*a+a2L-a2R);

        #b10
        bas[loc_idx(m,k,i+1,j),5+nbas*(k-1)]=-f*2;
        bas[loc_idx(m,umod(k-1,N),j,i+1),5+nbas*(k-1)]=-f*2;

        bas[loc_idx(m,k,i,j+1),5+nbas*(k-1)]=-f*a;

        bas[loc_idx(m,k,i-1,j),5+nbas*(k-1)]=2*f;
        bas[loc_idx(m,umod(k-1,N),j,i-1),5+nbas*(k-1)]=2*f;

        bas[loc_idx(m,umod(k-1,N),j+1,i),5+nbas*(k-1)]=-f*a;

        #b20
        bas[loc_idx(m,k,i+2,j),5+nbas*(k-1)]=-f;
        bas[loc_idx(m,umod(k-1,N),j,i+2),5+nbas*(k-1)]=-f;
        bas[loc_idx(m,k,i-2,j),5+nbas*(k-1)]=f+1;
        bas[loc_idx(m,umod(k-1,N),j,i-2),5+nbas*(k-1)]=f+1;

        g=1/(40*a+5*(a2L-a2R));
        #b11
        bas[loc_idx(m,k,i+1,j+1),5+nbas*(k-1)]=a*g*(3*a-2*a2R-10); #b11_1
        bas[loc_idx(m,k,i-1,j+1),5+nbas*(k-1)]=-a*g*(13*a-10-2*a2R); #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),5+nbas*(k-1)]=-a*g*(13*a-10-2*a2R); #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),5+nbas*(k-1)]=a*g*(3*a-2*a2R-10); #b11_0

        #b21_12

        bR=1/2*(g*a*(6*a2R-a3R-10));
        bL=1/2*(g*(16*a^2+(28*a2L+a3L-2*a2R+90)*a+4*(a2L+5/2)*(a2L-a2R)));

        bas[loc_idx(m,k,i+2,j+1),5+nbas*(k-1)]=bR; #b21_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+2),5+nbas*(k-1)]=bR;

        bas[loc_idx(m,k,i-2,j+1),5+nbas*(k-1)]=bL; #b21_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-2),5+nbas*(k-1)]=bL;

        #b31_13

        bbR=1/2*(g*3*a*a3R);
        bbL=1/2*(g*(-4*a2L^2+(2*a3L-32*a+4*a2R+10)*a2L+(13*a3L+80)*a-2*a2R*(a3L+5)));

        bas[loc_idx(m,k,i+3,j+1),5+nbas*(k-1)]=bbR; #b31_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+3),5+nbas*(k-1)]=bbR;

        bas[loc_idx(m,k,i-3,j+1),5+nbas*(k-1)]=bbL; #b31_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-3),5+nbas*(k-1)]=bbL;

        #B6

        #b10

        bas[loc_idx(m,k,i,j+1),6+nbas*(k-1)]=1;
        bas[loc_idx(m,umod(k-1,N),j+1,i),6+nbas*(k-1)]=-1;

        #b11
        bas[loc_idx(m,k,i+1,j+1),6+nbas*(k-1)]=1/2; #b11_1
        bas[loc_idx(m,k,i-1,j+1),6+nbas*(k-1)]=1/2; #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),6+nbas*(k-1)]=-1/2; #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),6+nbas*(k-1)]=-1/2; #b11_0 =#




    end

    for i in 1:size(bas,2)
       n=norm(bas[:,i]);
       bas[:,i]=bas[:,i]/2;
   end


    return bas


end


#=function g1basis_noncrossingRV_alternative(kn::Array, gluL::Array, gluR::Array, i::Int64, j::Int64, N::Int64)

    m, d = dim_deg(kn);
    bas=spzeros(N*m^2,4*N);

    a=gluR[1];
    a2R=gluR[2];
    a2L=-gluL[2];
    a3R=gluR[3];
    a3L=-gluL[1];

    #B1
    #b00=1
    for k in 1:N
        bas[loc_idx(m,k,i,j),1+4*(k-1)]=1;
        bas[loc_idx(m,umod(k-1,N),j,i),1+4*(k-1)]=1;

        #b10
        bas[loc_idx(m,k,i+1,j),1+4*(k-1)]=1;
        bas[loc_idx(m,umod(k-1,N),j,i+1),1+4*(k-1)]=1;

        bas[loc_idx(m,k,i,j+1),1+4*(k-1)]=1;

        bas[loc_idx(m,k,i-1,j),1+4*(k-1)]=1;
        bas[loc_idx(m,umod(k-1,N),j,i-1),1+4*(k-1)]=1;

        bas[loc_idx(m,umod(k-1,N),j+1,i),1+4*(k-1)]=1;

        #b20
        bas[loc_idx(m,k,i+2,j),1+4*(k-1)]=2/5;
        bas[loc_idx(m,umod(k-1,N),j,i+2),1+4*(k-1)]=2/5;
        bas[loc_idx(m,k,i-2,j),1+4*(k-1)]=2/5;
        bas[loc_idx(m,umod(k-1,N),j,i-2),1+4*(k-1)]=2/5;
        bas[loc_idx(m,k,i,j+2),1+4*(k-1)]=2/5;
        bas[loc_idx(m,umod(k-1,N),j+2,i),1+4*(k-1)]=2/5;

        #b11
        bas[loc_idx(m,k,i+1,j+1),1+4*(k-1)]=1-6/25*a; #b11_1
        bas[loc_idx(m,k,i-1,j+1),1+4*(k-1)]=1-6/25*a; #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),1+4*(k-1)]=1+6/25*a; #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),1+4*(k-1)]=1+6/25*a; #b11_0


        #b21_12
        bR=2/5*(1-3/5*a2R);
        bL=2/5*(1-3/5*a2L);
        bzero=2/5;

        bas[loc_idx(m,k,i+2,j+1),1+4*(k-1)]=bR; #b21_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+2),1+4*(k-1)]=bR;

        bas[loc_idx(m,k,i+1,j+2),1+4*(k-1)]=bzero; #b21_2
        bas[loc_idx(m,k,i-1,j+2),1+4*(k-1)]=bzero;

        bas[loc_idx(m,k,i-2,j+1),1+4*(k-1)]=bL; #b21_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-2),1+4*(k-1)]=bL;

        bas[loc_idx(m,umod(k-1,N),j+2,i+1),1+4*(k-1)]=bzero; #b21_0
        bas[loc_idx(m,umod(k-1,N),j+2,i-1),1+4*(k-1)]=bzero;

        #b31_13
        bbR=-6/25*a3R;
        bbL=-6/25*a3L;

        bas[loc_idx(m,k,i+3,j+1),1+4*(k-1)]=bbR; #b31_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+3),1+4*(k-1)]=bbR;

        bas[loc_idx(m,k,i-3,j+1),1+4*(k-1)]=bbL; #b31_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-3),1+4*(k-1)]=bbL;

        #B2
        #b00=0

        bas[loc_idx(m,k,i,j),2+4*(k-1)]=0;
        bas[loc_idx(m,umod(k-1,N),j,i),2+4*(k-1)]=0;

        #b10
        bas[loc_idx(m,k,i+1,j),2+4*(k-1)]=-1; #b10_1
        bas[loc_idx(m,umod(k-1,N),j,i+1),2+4*(k-1)]=-1;

        bas[loc_idx(m,k,i,j+1),2+4*(k-1)]=0; #b10_2

        bas[loc_idx(m,k,i-1,j),2+4*(k-1)]=1; #b10_3
        bas[loc_idx(m,umod(k-1,N),j,i-1),2+4*(k-1)]=1;

        bas[loc_idx(m,umod(k-1,N),j+1,i),2+4*(k-1)]=0; #b10_0

        #b20
        bas[loc_idx(m,k,i+2,j),2+4*(k-1)]=-1/2; #b20_1
        bas[loc_idx(m,umod(k-1,N),j,i+2),2+4*(k-1)]=-1/2;
        bas[loc_idx(m,k,i-2,j),2+4*(k-1)]=1/2; #b20_3
        bas[loc_idx(m,umod(k-1,N),j,i-2),2+4*(k-1)]=1/2;
        bas[loc_idx(m,k,i,j+2),2+4*(k-1)]=0; #b20_2
        bas[loc_idx(m,umod(k-1,N),j+2,i),2+4*(k-1)]=0; #b20_0

        #b11
        bas[loc_idx(m,k,i+1,j+1),2+4*(k-1)]=-1; #b11_1
        bas[loc_idx(m,k,i-1,j+1),2+4*(k-1)]=1; #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),2+4*(k-1)]=1; #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),2+4*(k-1)]=-1; #b11_0

        #b21_12
        bR=1/2*(-1+1/5*(3*a2R-1/2*a3R));
        bL=1/2*(1-1/5*(3*a2L-1/2*a3L));

        bas[loc_idx(m,k,i+2,j+1),2+4*(k-1)]=bR; #b21_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+2),2+4*(k-1)]=bR;

        bas[loc_idx(m,k,i-2,j+1),2+4*(k-1)]=bL; #b21_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-2),2+4*(k-1)]=bL;

        #b31_13
        bbR=3/10*a3R;
        bbL=-3/10*a3L;

        bas[loc_idx(m,k,i+3,j+1),2+4*(k-1)]=bbR; #b31_1
        bas[loc_idx(m,umod(k-1,N),j+1,i+3),2+4*(k-1)]=bbR;

        bas[loc_idx(m,k,i-3,j+1),2+4*(k-1)]=bbL; #b31_3
        bas[loc_idx(m,umod(k-1,N),j+1,i-3),2+4*(k-1)]=bbL;


        #B2
        #b00=0

        bas[loc_idx(m,k,i,j),3+4*(k-1)]=0;
        bas[loc_idx(m,umod(k-1,N),j,i),3+4*(k-1)]=0;

        #b10
        bas[loc_idx(m,k,i+1,j),3+4*(k-1)]=0; #b10_1
        bas[loc_idx(m,umod(k-1,N),j,i+1),3+4*(k-1)]=0;

        bas[loc_idx(m,k,i,j+1),3+4*(k-1)]=-1; #b10_2

        bas[loc_idx(m,k,i-1,j),3+4*(k-1)]=0; #b10_3
        bas[loc_idx(m,umod(k-1,N),j,i-1),3+4*(k-1)]=0;

        bas[loc_idx(m,umod(k-1,N),j+1,i),3+4*(k-1)]=1; #b10_0

        #b20
        bas[loc_idx(m,k,i+2,j),3+4*(k-1)]=0; #b20_1
        bas[loc_idx(m,umod(k-1,N),j,i+2),3+4*(k-1)]=0;
        bas[loc_idx(m,k,i-2,j),3+4*(k-1)]=0; #b20_3
        bas[loc_idx(m,umod(k-1,N),j,i-2),3+4*(k-1)]=0;
        bas[loc_idx(m,k,i,j+2),3+4*(k-1)]=-1/2; #b20_2
        bas[loc_idx(m,umod(k-1,N),j+2,i),3+4*(k-1)]=1/2; #b20_0

        #b11
        bas[loc_idx(m,k,i+1,j+1),3+4*(k-1)]=1; #b11_1
        bas[loc_idx(m,k,i-1,j+1),3+4*(k-1)]=-1; #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),3+4*(k-1)]=-1; #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),3+4*(k-1)]=1; #b11_0

        #b21_12
        bR=1/2*(-1+1/5*(3*a2R-1/2*a3R));
        bL=1/2*(1-1/5*(3*a2L-1/2*a3L));

        bas[loc_idx(m,k,i+1,j+2),3+4*(k-1)]=bR; #b21_2
        bas[loc_idx(m,k,i-1,j+2),3+4*(k-1)]=bR;

        bas[loc_idx(m,umod(k-1,N),j+2,i+1),3+4*(k-1)]=bL; #b21_4
        bas[loc_idx(m,umod(k-1,N),j+2,i-1),3+4*(k-1)]=bL;

        #b31_13
        bbR=3/10*a3R;
        bbL=-3/10*a3L;

        #bas[loc_idx(m,k,i+3,j+1),3+4*(k-1)]=bbR; #b31_1
        #bas[loc_idx(m,umod(k-1,N),j+1,i+3),3+4*(k-1)]=bbR;

        #bas[loc_idx(m,k,i-3,j+1),3+4*(k-1)]=bbL; #b31_3
        #bas[loc_idx(m,umod(k-1,N),j+1,i-3),3+4*(k-1)]=bbL;


        #B4
        #b11
        bas[loc_idx(m,k,i+1,j+1),4+4*(k-1)]=1; #b11_1
        bas[loc_idx(m,k,i-1,j+1),4+4*(k-1)]=-1; #b11_2
        bas[loc_idx(m,umod(k-1,N),j+1,i-1),4+4*(k-1)]=1; #b11_3
        bas[loc_idx(m,umod(k-1,N),j+1,i+1),4+4*(k-1)]=-1; #b11_0

    end


    return bas


end=#

function g1basis_noncrossingRV_alternative(kn::Array, gluL::Array, gluR::Array, i::Int64, j::Int64, N::Int64)

    m, d = dim_deg(kn);
    bas=spzeros(N*m^2,4*N);

    a=gluR[1];
    a2R=gluR[2];
    a2L=-gluL[2];
    a3R=gluR[3];
    a3L=-gluL[1];

    #B1
    #b00=1
    for k in 1:N
        bas[loc_idx(m,k,i,j),1+4*(k-1)]=1/4;
        bas[loc_idx(m,umod(k-1,N),j,i),1+4*(k-1)]=1/4;

        #b10
        bas[loc_idx(m,k,i+1,j),1+4*(k-1)]=1/2;
        bas[loc_idx(m,umod(k-1,N),j,i+1),1+4*(k-1)]=1/2;

        bas[loc_idx(m,k,i,j+1),1+4*(k-1)]=1/2;


        #b11
        bas[loc_idx(m,k,i+1,j+1),1+4*(k-1)]=1; #b11_1



        #B2
        #b00=0

        bas[loc_idx(m,k,i,j),2+4*(k-1)]=1/4;
        bas[loc_idx(m,umod(k-1,N),j,i),2+4*(k-1)]=1/4;

        #b10

        bas[loc_idx(m,k,i,j+1),2+4*(k-1)]=1/2; #b10_2

        bas[loc_idx(m,k,i-1,j),2+4*(k-1)]=1/2; #b10_3
        bas[loc_idx(m,umod(k-1,N),j,i-1),2+4*(k-1)]=1/2;


        #b11

        bas[loc_idx(m,k,i-1,j+1),2+4*(k-1)]=1; #b11_2

        #B3
        #b00=0

        bas[loc_idx(m,k,i,j),3+4*(k-1)]=1/4;
        bas[loc_idx(m,umod(k-1,N),j,i),3+4*(k-1)]=1/4;

        #b10


        bas[loc_idx(m,k,i-1,j),3+4*(k-1)]=1/2; #b10_3
        bas[loc_idx(m,umod(k-1,N),j,i-1),3+4*(k-1)]=1/2;

        bas[loc_idx(m,umod(k-1,N),j+1,i),3+4*(k-1)]=1/2; #b10_0



        #b11

        bas[loc_idx(m,umod(k-1,N),j+1,i-1),3+4*(k-1)]=1; #b11_3


        #B4
        #b00=0

        bas[loc_idx(m,k,i,j),4+4*(k-1)]=1/4;
        bas[loc_idx(m,umod(k-1,N),j,i),4+4*(k-1)]=1/4;

        #b10


        bas[loc_idx(m,k,i+1,j),4+4*(k-1)]=1/2; #b10_3
        bas[loc_idx(m,umod(k-1,N),j,i+1),4+4*(k-1)]=1/2;

        bas[loc_idx(m,umod(k-1,N),j+1,i),4+4*(k-1)]=1/2; #b10_0



        #b11

        bas[loc_idx(m,umod(k-1,N),j+1,i+1),4+4*(k-1)]=1; #b11_3
    end


    return bas


end


function g1basis_EVedge_gluing(kn::Array,nbsubd::Int64)

    m, d = dim_deg(kn)
    bas=spzeros(2*m^2,2^(nbsubd+1))

    if nbsubd==0
        j=1
        for i in 3:m-2
            bas[loc_idx(m,1,i,2),j]=1; #b21
            bas[loc_idx(m,2,2,i),j]=-1;
            j+=1
        end
    else

        m0=6;
        excl=[];

        for k in 1:2^nbsubd -1
            p=m0+(m0-1)*(k-1);
            push!(excl,p-1);
            push!(excl,p);
            push!(excl,p+1);
        end

        j=1

        for i in 3:m-2
            if !(i in excl)
                bas[loc_idx(m,1,i,2),j]=1; #b21
                bas[loc_idx(m,2,2,i),j]=-1;
                j+=1
            end
        end
    end

    return bas;
end


function g1basis_Redge_gluing(kn::Array,nbsubd::Int64)

    m, d = dim_deg(kn)
    bas=spzeros(2*m^2,2^(nbsubd+2))

    if nbsubd==0
        j=1
        for i in 3:m-2
            bas[loc_idx(m,1,i,1),j]=1/2; #b21
            bas[loc_idx(m,2,1,i),j]=1/2;
            bas[loc_idx(m,1,i,2),j]=1;
            j+=1
            bas[loc_idx(m,1,i,1),j]=1/2; #b21
            bas[loc_idx(m,2,1,i),j]=1/2;
            bas[loc_idx(m,2,2,i),j]=1;
            j+=1
        end
    else

        m0=6;
        excl=[];

        for k in 1:2^nbsubd -1
            p=m0+(m0-1)*(k-1);
            push!(excl,p-1);
            push!(excl,p);
            push!(excl,p+1);
        end

        j=1

        for i in 3:m-2
            if !(i in excl)
                bas[loc_idx(m,1,i,1),j]=1/2; #b21
                bas[loc_idx(m,2,1,i),j]=1/2;
                bas[loc_idx(m,1,i,2),j]=1;
                j+=1
                bas[loc_idx(m,1,i,1),j]=1/2; #b21
                bas[loc_idx(m,2,1,i),j]=1/2;
                bas[loc_idx(m,2,2,i),j]=1;
                j+=1
            end
        end
    end

    return bas;
end
