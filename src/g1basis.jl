export g1basis
"""
This function takes in input a quad mesh and returns a sparse matrix containing the coefficients defining a set of G1 biquintic basis functions on the input mesh.

 - ncols of the sparse matrix gives the dimension of the spline space,
 - nrows is the total number of control points in the mesh i.e. nfaces*36
"""
function g1basis(m::HMesh)
    basis=g1basis_bezier(m);
    return basis
end


function g1basis(m::Mesh)
    hm=hmesh(m)
    g1basis(hm);
end


function g1basis_EV(N::Int64)

#dim=12;

bas=spzeros(36*N,N+3);
#I fix b00=1 and all the other vertices to 0
a=2*cos(2*pi/N);

#B1

#b00

VAL=0.5;

for i in 0:N-1
    bas[i*36+1,1]=VAL; #b00=1
end

#b10

c=zeros(N); #Right hand of the system
c[1]=1; #Vector for the circulant matrix
c[N-1]=1;
c[N]=-a;
T=Circulant(c);
T=convert(Matrix,T);
k=nullspace(T);
F=[T;k[:,1]';k[:,2]'];
RHS=zeros(N+2);
RHS[1:N].=(2-a)*VAL;
#RHS[N+1]=-1.5;
#RHS[N+2]=-1.5;
sol=F\RHS;
push!(sol,sol[1]);
for i in 0:N-1
    bas[i*36+2,1]=sol[i+1]; #b10
    bas[i*36+7,1]=sol[i+2]; #b01
end

#b20

for i in 0:N-1
    bas[i*36+3,1]=1/2*bas[i*36+2,1]-VAL/10; #b20
    bas[i*36+13,1]=1/2*bas[i*36+7,1]-VAL/10; #b02
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
    for i in 0:N-1
        c[i+1]=1/5*(a*VAL-5*(a-2)*bas[i*36+2,1]+4*a*bas[i*36+3,1]);
    end
    ff=[c;0];
    sol=Z\ff;

    for i in 0:N-1
        bas[i*36+8,1]=sol[i+1]; #b11
    end
else
    for i in 0:N-1
        c[i+1]=1/5*(a*VAL-5*(a-2)*bas[i*36+2,1]+4*a*bas[i*36+3,1]);
    end
    sol=T\c;

    for i in 0:N-1
        bas[i*36+8,1]=sol[i+1]; #b11
    end
end


#b21-b12 and b31-b13

for i in 0:N-1
    R=1/10*(-a*VAL+5*a*bas[i*36+2,1]-10*(a-2)*bas[i*36+3,1]);
    bas[i%N*36+9,1]=1/2*R; #b21
    bas[(N-1+i)%N*36+14,1]=1/2*R; #b12
    R1=1/10*(a*VAL-5*a*bas[i*36+2,1]+10*a*bas[i*36+3,1]);
    bas[i%N*36+10,1]=1/2*R1; #b31
    bas[(N-1+i)%N*36+20,1]=1/2*R1; #b13
end


#B2-B3

#b00

for i in 0:N-1
    bas[i*36+1,2]=0; #b00=0
    bas[i*36+1,3]=0; #b00=0
end

#b10

v=zeros(N);
v[1]=1;
v[N-1]=1;
v[N]=-a;

H=Circulant(v);
H=convert(Matrix,H);
sol=nullspace(H);
sol1=sol[:,1];
sol1=sol1./maximum(sol1);
push!(sol1,sol1[1]);
sol2=sol[:,2];
sol2=sol2./maximum(sol2);
push!(sol2,sol2[1]);

for i in 0:N-1
    bas[i*36+2,2]=sol1[i+1]; #b10
    bas[i*36+7,2]=sol1[i+2]; #b01
    bas[i*36+2,3]=sol2[i+1]; #b10
    bas[i*36+7,3]=sol2[i+2]; #b01
end

#b20

for i in 0:N-1
    bas[i*36+3,2]=1/2*bas[i*36+2,2]; #b20
    bas[i*36+13,2]=1/2*bas[i*36+7,2]; #b02
    bas[i*36+3,3]=1/2*bas[i*36+2,3]; #b20
    bas[i*36+13,3]=1/2*bas[i*36+7,3]; #b02
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
    for i in 0:N-1
        c1[i+1]=(1/5)*(-5*(a-2)*bas[i*36+2,2]+4*a*bas[i*36+3,2]);
        c2[i+1]=(1/5)*(-5*(a-2)*bas[i*36+2,3]+4*a*bas[i*36+3,3]);
    end
    cc1=[c1;0];
    cc2=[c2;0];
    sol1_11=Z\cc1;
    sol2_11=Z\cc2;

    for i in 0:N-1
        bas[i*36+8,2]=sol1_11[i+1];
        bas[i*36+8,3]=sol2_11[i+1]; #b11
    end

else

    for i in 0:N-1
        c1[i+1]=(1/5)*(-5*(a-2)*bas[i*36+2,2]+4*a*bas[i*36+3,2]);
        c2[i+1]=(1/5)*(-5*(a-2)*bas[i*36+2,3]+4*a*bas[i*36+3,3]);
    end
    sol1_11=T\c1;
    sol2_11=T\c2;

    for i in 0:N-1
        bas[i*36+8,2]=sol1_11[i+1];
        bas[i*36+8,3]=sol2_11[i+1]; #b11
    end
end


#b21-b12 and b31-b13


for i in 0:N-1
    R1=(1/10)*(5*a*bas[i*36+2,2]-10*(a-2)*bas[i*36+3,2]);
    bas[i%N*36+9,2]=1/2*R1; #b21
    bas[(N-1+i)%N*36+14,2]=1/2*R1; #b12
    R2=(1/10)*(5*a*bas[i*36+2,3]-10*(a-2)*bas[i*36+3,3]);
    bas[i%N*36+9,3]=1/2*R2; #b21
    bas[(N-1+i)%N*36+14,3]=1/2*R2; #b12

    R3=(1/10)*(-5*a*bas[i*36+2,2]+10*a*bas[i*36+3,2]);
    bas[i%N*36+10,2]=1/2*R3; #b31
    bas[(N-1+i)%N*36+20,2]=1/2*R3; #b13
    R4=(1/10)*(-5*a*bas[i*36+2,3]+10*a*bas[i*36+3,3]);
    bas[i%N*36+10,3]=1/2*R4; #b31
    bas[(N-1+i)%N*36+20,3]=1/2*R4; #b13
end



#B4-B6

#b11

for i in 0:N-1

    bas[i*36+8,4+i]=VAL; #b11^0=1

    bas[i%N*36+3,4+i]=(5/(4*a))*bas[i*36+8,4+i]*VAL; #b20
    bas[(N-1+i)%N*36+13,4+i]=(5/(4*a))*bas[i*36+8,4+i]*VAL; #b02
    bas[i%N*36+13,4+i]=(5/(4*a))*bas[i*36+8,4+i]*VAL; #b20
    bas[(i+1)%N*36+3,4+i]=(5/(4*a))*bas[i*36+8,4+i]*VAL; #b02

#b30

    bas[i%N*36+4,4+i]=bas[i*36+3,4+i];
    bas[(N-1+i)%N*36+19,4+i]=bas[i*36+3,4+i];
    bas[i%N*36+19,4+i]=bas[i*36+3,4+i];
    bas[(i+1)%N*36+4,4+i]=bas[i*36+3,4+i];

#b21-b12 b31-b13

for k in 0:1
    R=-(a-2)*bas[(i+k)%N*36+3,4+i]+(3/5)*a*bas[i*36+4,4+i];
    bas[(i+k)%N*36+9,4+i]=1/2*R; #b21
    bas[(N-1+k+i)%N*36+14,4+i]=1/2*R; #b12
    R=a*bas[(i+k)%N*36+3,4+i]-(a-2)*bas[i*36+4,4+i];
    bas[(i+k)%N*36+10,4+i]=1/2*R; #b31
    bas[(N-1+k+i)%N*36+20,4+i]=1/2*R; #b13
end

#=
for i in 0:N-1
        bas[i%N*36+9,7+i]=1; #b21
        bas[(N-1+i)%N*36+14,7+i]=-1; #b12
end

for i in 0:N-1
        bas[i%N*36+10,10+i]=1; #b31
        bas[(N-1+i)%N*36+20,10+i]=-1; #b13
end

=#

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

function g1basis_EVedge()

    bas=spzeros(2*36,1)

    #B7-B9


    bas[9,1]=1; #b21
    bas[36+14,1]=-1; #b12


    #B10-B12
    #=
    for i in 0:N-1
            bas[i%N*36+10,10+i]=1; #b21
            bas[(N-1+i)%N*36+20,10+i]=-1; #b12
    end
    =#
    return bas
end

function g1basis_RV(v::Array)

#p is used to modify even b30 if the edge is attached to an EV

bas=spzeros(4*36,4)
sol10=[1/2,1/2,0,0];
sol11=[1,0,0,0];

for k in 1:4
    for i in 0:3
        bas[i*36+1,k]=1/4; #b00
    end
    for i in 0:3
        bas[i*36+2,k]=sol10[i%4+1]; #b10
        bas[i*36+7,k]=sol10[(i+1)%4+1]; #b01
        bas[i*36+8,k]=sol11[i+1]; #b11
    end
    sol10=vectorshift(sol10,1);
    sol11=vectorshift(sol11,1);
end

for h in 0:3
    if v[1]!=0
        N=v[1];
        a=2*cos(2*pi/N);
        bas[(h%4)*36+3,h+1]=9/40;
        bas[(h+3)%4*36+13,h+1]=bas[(h%4)*36+3,h+1];
        bas[(h%4)*36+9,h+1]=9/40-a/80;
        bas[(h+3)%4*36+14,h+1]=bas[(h%4)*36+9,h+1];
        bas[(h%4)*36+10,h+1]=27/400*a;
        bas[(h+3)%4*36+20,h+1]=bas[(h%4)*36+10,h+1];
    end
    if v[2]!=0
        N2=v[2];
        a=2*cos(2*pi/N2);
        bas[(h%4)*36+13,h+1]=9/40;
        bas[(h+1)%4*36+3,h+1]=bas[(h%4)*36+13,h+1];
        bas[(h%4)*36+14,h+1]=9/40-a/80;
        bas[(h+1)%4*36+9,h+1]=bas[(h%4)*36+14,h+1];
        bas[(h%4)*36+20,h+1]=27/400*a;
        bas[(h+1)%4*36+10,h+1]=bas[(h%4)*36+20,h+1];
    end
    if v[3]!=0
        N3=v[3];
        a=2*cos(2*pi/N3);
        bas[(h+1)%4*36+13,h+1]=-1/40;
        bas[(h+2)%4*36+3,h+1]=bas[(h+1)%4*36+13,h+1];
        bas[(h+1)%4*36+14,h+1]=-(2-a)/80;
        bas[(h+2)%4*36+9,h+1]=bas[(h+1)%4*36+14,h+1];
        bas[(h+1)%4*36+20,h+1]=-3*a/400;
        bas[(h+2)%4*36+10,h+1]=bas[(h+1)%4*36+20,h+1];
    end
    if v[4]!=0
        N4=v[4];
        a=2*cos(2*pi/N4);
        bas[(h+2)%4*36+13,h+1]=-1/40;
        bas[(h+3)%4*36+3,h+1]=bas[(h+2)%4*36+13,h+1];
        bas[(h+2)%4*36+14,h+1]=-(2-a)/80;
        bas[(h+3)%4*36+9,h+1]=bas[(h+2)%4*36+14,h+1];
        bas[(h+2)%4*36+20,h+1]=-3*a/400;
        bas[(h+3)%4*36+10,h+1]=bas[(h+2)%4*36+20,h+1];
    end
    v=vectorshift(v,-1)
end

#=
if p!=0
    bas[3,p]=1/2*bas[2,p]-1/10*bas[1,p];
    bas[4*36+13,p]=bas[3,p];
end
=#

return bas
end

function g1basis_Redge()

bas=spzeros(2*36,2)

bas[3,1]=1/2;
bas[9,1]=1; #b21
bas[36+13,1]=1/2;

bas[3,2]=1/2;
bas[36+13,2]=1/2; #b21
bas[36+14,2]=1;

return bas

end

function g1basis_borderRV(N::Int64=0)

bas=spzeros(2*36,4)

bas[1,1]=1/2
bas[36+1,1]=1/2
bas[2,1]=1

bas[7,2]=1/2
bas[36+2,2]=1/2
bas[8,2]=1

bas[7,3]=1/2
bas[36+2,3]=1/2
bas[36+8,3]=1

bas[1,4]=1/2
bas[36+1,4]=1/2
bas[36+7,4]=1

if N!=0
    a=2*cos(2*pi/N);
    for k in 1:4
        bas[13,k]=1/2*bas[7,k]-1/10*bas[1,k];
        bas[36+3,k]=bas[13,k];
        bas[14,k]=1/2*(2-a)*bas[13,k]+1/5*a*bas[7,k];
        bas[36+9,k]=bas[14,k];
        bas[20,k]=3/10*a*bas[13,k];
        bas[36+10,k]=bas[20,k];
    end
end

return bas

end

function g1basis_Faces()

bas=spzeros(1*36,4)

bas[15,1]=1;
bas[16,2]=1;
bas[21,3]=1;
bas[22,4]=1;

return bas

end

function g1basis_Corners()

bas=spzeros(1*36,4)

bas[1,1]=1;
bas[2,2]=1;
bas[7,3]=1;
bas[8,4]=1;

return bas

end

function g1basis_RBedge()

bas=spzeros(1*36,4)

bas[3,1]=1;
bas[4,2]=1;
bas[9,3]=1;
bas[10,4]=1;

return bas

end
