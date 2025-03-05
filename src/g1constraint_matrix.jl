#=
This function creates the smoothing matrix to create a G1 smooth patch around extraordinary vertices of any valence
 - N--> Valence of the EV
 - S--> Solving strategy NCS-S (non circulant system+symm), NCS-AS(non circulant system+asymm), CS-AS(circulant system+asymm), CS-S(circulant system+symm)
T-->Type of mesh face
The output matrix B has smoothing masks ordered as follow

  Control points labelling for inner EVs
      Biquintic patch
```
    28 31 34 | 21 20 19
    29 32 35 | 24 23 22
    30 33 36 | 25 26 27
    -------------------       B=[1|2|3|...|34|35|36] columns OF B
    7  8  9  | 18 15 12
    4  5  6  | 17 14 11
    1  2  3  | 16 13 10
```
 For regular vertices N=4 it will return ACC3 smoothing masks
=#
function G1matrix(N::Int, S::String="CS-S", T::String="INNER")


if S=="NCS-S"
    circulant_system="no";
    symm="yes";
elseif S=="NCS-AS"
    circulant_system="no";
    symm="no";
elseif S=="CS-AS"
    circulant_system="yes";
    symm="no";
elseif S=="CS-S"
    circulant_system="yes";
    symm="yes";
end

if T=="INNER"
    if N==4
        B=acc3bigmatrix(N,"INNEREV");
    else
        M2=degree_elevate_mask(N,"INNEREV",0); #1
        B=ordermatrix(M2,N,"INNER");
        a=2*cos(2*pi/N);

        if N!=3 #Compute the new M10

            c=zeros(N); #Right hand of the system
            c[1]=-a; #Vector for the circulant matrix
            c[2]=1;
            c[N]=1;
            T=Circulant(c);
            T=convert(Matrix,T);
            k=nullspace(T);
            F=[T;k[:,1]';k[:,2]'];
            b1=dot(M2[1,2][2:N+1],k[:,1]);
            b2=dot(M2[1,2][2:N+1],k[:,2]);
            w=M2[1,1][1];
            c=(2-a)*M2[1,1][2:N+1];
            c=[c;b1;b2];
            x=F\c;
            b1=dot(M2[1,2][N+2:2*N+1],k[:,1]);
            b2=dot(M2[1,2][N+2:2*N+1],k[:,2]);
            c=(2-a)*M2[1,1][N+2:2*N+1];
            c=[c;b1;b2];
            y=F\c;
            M10new=[w;x;y];

            M10new=normalizemask(M10new);
            B[1:2*N+1,2]=M10new;

            M01=completeshift(M10new,N,1,"INNER");
            B[:,4]=M01;

            #Compute new M20
            M20new=(1/10)*(M2[1,6]-M2[1,1]+5*B[:,2]+10*M2[1,4]-5*M2[1,5]);

            M20new=normalizemask(M20new);
            B[:,3]=M20new;

            M02=completeshift(M20new,N,1,"INNER");
            B[:,7]=M02;
        else
            M20new=(1/10)*(M2[1,6]-M2[1,1]+5*M2[1,2]+10*M2[1,4]-5*M2[1,5]);

            M20new=normalizemask(M20new);
            B[:,3]=M20new;

            M02=completeshift(M20new,N,1,"INNER");
            B[:,7]=M02;
        end

        #M11

        if N%2==0   #Even case I modify here M20VER-->M11-->M30
            if circulant_system=="yes"
                v=zeros(N); #Circulant matrix vector
                k=zeros(1,N); #It will be the kernel vector -1 1 -1 ... 1
                for i in 1:N
                    k[i]=(-1)^i;
                end
                v[1]=1; #Vector for the circulant matrix
                v[N]=1;
                A=convert(Matrix,Circulant(v)); #Circulant matrix for the M11 LS system
                Z=[A;k];

                #M20VER

                f=zeros(N);
                for i in 1:N
                    f[i]=dot(B[2:N+1,3],A[i,:]);
                end
                ff=[f;0];
                xM20new=Z\ff;
                for i in 1:N
                    f[i]=dot(B[N+2:2*N+1,3],A[i,:]);
                end
                ff=[f;0];
                yM20new=Z\ff;
                w=B[1,3];
                M20even_new=[w;xM20new;yM20new];

                M20even_new=normalizemask(M20even_new);
                B[1:2*N+1,3]=M20even_new;

                M02=completeshift(M20even_new,N,1,"INNER");
                B[:,7]=M02;

                #M11

                w=(1/10)*(a*M2[1,1][1]+(10-5*a)*B[1,2]+4*a*M20even_new[1]);
                c=(1/5)*(a*M2[1,1][2:N+1]+(10-5*a)*B[2:N+1,2]+4*a*B[2:N+1,3]);
                c=[c;dot(k,M2[2,2][2:N+1])];
                x=Z\c; #First neighborhood solution
                cc=(1/5)*(a*M2[1,1][N+2:2*N+1]+(10-5*a)*B[N+2:2*N+1,2]+4*a*B[N+2:2*N+1,3]);
                cc=[cc;dot(k,M2[2,2][N+2:2*N+1])];
                y=Z\cc; #Second neighborhood solution
                M11new=[w;x;y];

                M11new=normalizemask(M11new);
                B[1:2*N+1,5]=M11new;

                #M30

                M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);
                M30new=normalizemask(M30new);

                B[:,16]=M30new;

                M03new=completeshift(M30new,N,1,"INNER");

                B[:,30]=M03new;

            else

                #Compute M20 from M11 and then M30

                #M20

                M11_0=completeshift(M2[2,2],N,-1,"INNER");
                M20new=(1/(4*a))*(5*M2[2,2]-a*M2[1,1]+5*(a-2)*B[:,2]+5*M11_0);

                M20new=normalizemask(M20new);
                B[:,3]=M20new;

                M02=completeshift(M20new,N,1,"INNER");
                B[:,7]=M02;

                #M30

                M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);

                M30new=normalizemask(M30new);
                B[:,16]=M30new;

                M03new=completeshift(M30new,N,1,"INNER");

                M03new=normalizemask(M03new);
                B[:,30]=M03new;

            end

            elseif N%2!=0 #Odd case
                if circulant_system=="yes"

                    #M11

                    v=zeros(N);
                    v[1]=1;
                    v[N]=1;

                    #Quadratic glueing data
                    w=(1/10)*(a*M2[1,1][1]+(10-5*a)*B[1,2]+4*a*B[1,3]);
                    A=convert(Matrix,Circulant(v));
                    c=(1/5)*(a*M2[1,1][2:N+1]+(10-5*a)*B[2:N+1,2]+4*a*B[2:N+1,3]);
                    x=A\c;
                    c=(1/5)*(a*M2[1,1][N+2:2*N+1]+(10-5*a)*B[N+2:2*N+1,2]+4*a*B[N+2:2*N+1,3]);
                    y=A\c;
                    M11new=[w;x;y];

                    M11new=normalizemask(M11new);
                    B[1:2*N+1,5]=M11new;

                else

                    #M20

                    M11_0=completeshift(M2[2,2],N,-1,"INNER");
                    M20new=(1/(4*a))*(5*M2[2,2]-a*M2[1,1]+5*(a-2)*B[:,2]+5*M11_0);

                    M20new=normalizemask(M20new);
                    B[:,3]=M20new;

                    M02=completeshift(M20new,N,1,"INNER");
                    B[:,7]=M02;

                    #M30

                    M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);

                    M30new=normalizemask(M30new);
                    B[:,16]=M30new;

                    M03new=completeshift(M30new,N,1,"INNER");
                    B[:,30]=M03new;

                end
            end

            #M21 non symmetric

            if symm=="no"

                M12_0=completeshift(M2[3,2],N,-1,"INNER");
                M21=(1/10)*(-a*M2[1,1]+5*a*B[:,2]+(20-10*a)*B[:,3]+6*a*B[:,16]-10*M12_0);

                M21=normalizemask(M21);
                B[:,6]=M21;

            #M21 symmetric

            else

                R=(1/10)*(-a*M2[1,1]+5*a*B[:,2]+(20-10*a)*B[:,3]+6*a*B[:,16]);
                A12_0=completeshift(M2[3,2],N,-1,"INNER");
                A21=M2[2,3];
                theta=(A21-A12_0);
                M21=1/2*(R+theta);

                M21=normalizemask(M21);
                B[:,6]=M21;

                M12_0=1/2*(R-theta);

                M12_0=normalizemask(M12_0);

                M12=completeshift(M12_0,N,1,"INNER");
                B[:,8]=M12;

            end

        #M31 non symmetric

        if symm=="no"

            M13_0=completeshift(M2[4,2],N,-1,"INNER");
            M31=(1/10)*(a*M2[1,1]-5*a*B[:,2]+10*a*B[:,3]+(20-10*a)*B[:,16]+4*a*M2[1,5]-10*M13_0);

            M31=normalizemask(M31);
            B[:,17]=M31;

        #M31 symmetric

        else

            R=(1/10)*(a*M2[1,1]-5*a*B[:,2]+10*a*B[:,3]+(20-10*a)*B[:,16]+4*a*M2[1,5]);
            A13_0=completeshift(M2[4,2],N,-1,"INNER");
            A31=M2[2,4];
            theta=(A31-A13_0);
            M31=1/2*(R+theta);

            M31=normalizemask(M31);
            B[:,17]=M31;

            M13_0=1/2*(R-theta);
            M13_0=normalizemask(M13_0);
            M13=completeshift(M13_0,N,1,"INNER");
            B[:,33]=M13;

        end
    end

elseif T=="EVREGBORDER"

    M2=degree_elevate_mask(N, "CORNER",0); #6
    B=ordermatrix(M2,N,"EVREGBORDER");
    a=2*cos(2*pi/N);

    if N!=3 #Compute the new M10

        c=zeros(N); #Right hand of the system
        c[1]=-a; #Vector for the circulant matrix
        c[2]=1;
        c[N]=1;
        T=Circulant(c);
        T=convert(Matrix,T);
        k=nullspace(T);
        F=[T;k[:,1]';k[:,2]'];
        b1=dot(M2[1,2][2:N+1],k[:,1]);
        b2=dot(M2[1,2][2:N+1],k[:,2]);
        w=M2[1,1][1];
        c=(2-a)*M2[1,1][2:N+1];
        c=[c;b1;b2];
        x=F\c;
        b1=dot(M2[1,2][N+2:2*N+1],k[:,1]);
        b2=dot(M2[1,2][N+2:2*N+1],k[:,2]);
        c=(2-a)*M2[1,1][N+2:2*N+1];
        c=[c;b1;b2];
        y=F\c;
        M10new=[w;x;y];

        M10new=normalizemask(M10new);
        B[1:2*N+1,2]=M10new;

        M01=completeshift(M10new,N,1,"EVREGBORDER");
        B[:,4]=M01;

        #Compute new M20
        M20new=(1/10)*(M2[1,6]-M2[1,1]+5*B[:,2]+10*M2[1,4]-5*M2[1,5]);

        M20new=normalizemask(M20new);
        B[:,3]=M20new;

        M02=completeshift(M20new,N,1,"EVREGBORDER");
        B[:,7]=M02;
    else
        M20new=(1/10)*(M2[1,6]-M2[1,1]+5*M2[1,2]+10*M2[1,4]-5*M2[1,5]);

        M20new=normalizemask(M20new);
        B[:,3]=M20new;

        M02=completeshift(M20new,N,1,"EVREGBORDER");
        B[:,7]=M02;
    end

    #M11

    if N%2==0   #Even case I modify here M20VER-->M11-->M30
        if circulant_system=="yes"
            v=zeros(N); #Circulant matrix vector
            k=zeros(1,N); #It will be the kernel vector -1 1 -1 ... 1
            for i in 1:N
                k[i]=(-1)^i;
            end
            v[1]=1; #Vector for the circulant matrix
            v[N]=1;
            A=convert(Matrix,Circulant(v)); #Circulant matrix for the M11 LS system
            Z=[A;k];

            #M20VER

            f=zeros(N);
            for i in 1:N
                f[i]=dot(B[2:N+1,3],A[i,:]);
            end
            ff=[f;0];
            xM20new=Z\ff;
            for i in 1:N
                f[i]=dot(B[N+2:2*N+1,3],A[i,:]);
            end
            ff=[f;0];
            yM20new=Z\ff;
            w=B[1,3];
            M20even_new=[w;xM20new;yM20new];

            M20even_new=normalizemask(M20even_new);
            B[1:2*N+1,3]=M20even_new;

            M02=completeshift(M20even_new,N,1,"EVREGBORDER");
            B[:,7]=M02;

            #M11

            w=(1/10)*(a*M2[1,1][1]+(10-5*a)*B[1,2]+4*a*M20even_new[1]);
            c=(1/5)*(a*M2[1,1][2:N+1]+(10-5*a)*B[2:N+1,2]+4*a*B[2:N+1,3]);
            c=[c;dot(k,M2[2,2][2:N+1])];
            x=Z\c; #First neighborhood solution
            cc=(1/5)*(a*M2[1,1][N+2:2*N+1]+(10-5*a)*B[N+2:2*N+1,2]+4*a*B[N+2:2*N+1,3]);
            cc=[cc;dot(k,M2[2,2][N+2:2*N+1])];
            y=Z\cc; #Second neighborhood solution
            M11new=[w;x;y];

            M11new=normalizemask(M11new);
            B[1:2*N+1,5]=M11new;

            #M30

            M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);
            M30new=normalizemask(M30new);

            B[:,16]=M30new;

            M03new=completeshift(M30new,N,1,"EVREGBORDER");

            B[:,30]=M03new;

        else

            #Compute M20 from M11 and then M30

            #M20

            M11_0=completeshift(M2[2,2],N,-1,"EVREGBORDER");
            M20new=(1/(4*a))*(5*M2[2,2]-a*M2[1,1]+5*(a-2)*B[:,2]+5*M11_0);

            M20new=normalizemask(M20new);
            B[:,3]=M20new;

            M02=completeshift(M20new,N,1,"EVREGBORDER");
            B[:,7]=M02;

            #M30

            M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);

            M30new=normalizemask(M30new);
            B[:,16]=M30new;

            M03new=completeshift(M30new,N,1,"EVREGBORDER");

            M03new=normalizemask(M03new);
            B[:,30]=M03new;

        end

        elseif N%2!=0 #Odd case
            if circulant_system=="yes"

                #M11

                v=zeros(N);
                v[1]=1;
                v[N]=1;

                #Quadratic glueing data
                w=(1/10)*(a*M2[1,1][1]+(10-5*a)*B[1,2]+4*a*B[1,3]);
                A=convert(Matrix,Circulant(v));
                c=(1/5)*(a*M2[1,1][2:N+1]+(10-5*a)*B[2:N+1,2]+4*a*B[2:N+1,3]);
                x=A\c;
                c=(1/5)*(a*M2[1,1][N+2:2*N+1]+(10-5*a)*B[N+2:2*N+1,2]+4*a*B[N+2:2*N+1,3]);
                y=A\c;
                M11new=[w;x;y];

                M11new=normalizemask(M11new);
                B[1:2*N+1,5]=M11new;

            else

                #M20

                M11_0=completeshift(M2[2,2],N,-1,"EVREGBORDER");
                M20new=(1/(4*a))*(5*M2[2,2]-a*M2[1,1]+5*(a-2)*B[:,2]+5*M11_0);

                M20new=normalizemask(M20new);
                B[:,3]=M20new;

                M02=completeshift(M20new,N,1,"EVREGBORDER");
                B[:,7]=M02;

                #M30

                M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);

                M30new=normalizemask(M30new);
                B[:,16]=M30new;

                M03new=completeshift(M30new,N,1,"EVREGBORDER");
                B[:,30]=M03new;

            end
        end

        #M21 non symmetric

        if symm=="no"

            M12_0=completeshift(M2[3,2],N,-1,"EVREGBORDER");
            M21=(1/10)*(-a*M2[1,1]+5*a*B[:,2]+(20-10*a)*B[:,3]+6*a*B[:,16]-10*M12_0);

            M21=normalizemask(M21);
            B[:,6]=M21;

        #M21 symmetric

        else

            R=(1/10)*(-a*M2[1,1]+5*a*B[:,2]+(20-10*a)*B[:,3]+6*a*B[:,16]);
            A12_0=completeshift(M2[3,2],N,-1,"EVREGBORDER");
            A21=M2[2,3];
            theta=(A21-A12_0);
            M21=1/2*(R+theta);

            M21=normalizemask(M21);
            B[:,6]=M21;

            M12_0=1/2*(R-theta);

            M12_0=normalizemask(M12_0);

            M12=completeshift(M12_0,N,1,"EVREGBORDER");
            B[:,8]=M12;

        end

    #M31 non symmetric

    if symm=="no"

        M13_0=completeshift(M2[4,2],N,-1,"EVREGBORDER");
        M31=(1/10)*(a*M2[1,1]-5*a*B[:,2]+10*a*B[:,3]+(20-10*a)*B[:,16]+4*a*M2[1,5]-10*M13_0);

        M31=normalizemask(M31);
        B[:,17]=M31;

    #M31 symmetric

    else

        R=(1/10)*(a*M2[1,1]-5*a*B[:,2]+10*a*B[:,3]+(20-10*a)*B[:,16]+4*a*M2[1,5]);
        A13_0=completeshift(M2[4,2],N,-1,"EVREGBORDER");
        A31=M2[2,4];
        theta=(A31-A13_0);
        M31=1/2*(R+theta);

        M31=normalizemask(M31);
        B[:,17]=M31;

        M13_0=1/2*(R-theta);
        M13_0=normalizemask(M13_0);
        M13=completeshift(M13_0,N,1,"EVREGBORDER");
        B[:,33]=M13;

    end

elseif T=="EVREGBORDER1"

    M2=degree_elevate_mask(N,"BORDER1",0);#2
    B=ordermatrix(M2,N,"EVREGBORDER1");
    a=2*cos(2*pi/N);

    if N!=3 #Compute the new M10

        c=zeros(N); #Right hand of the system
        c[1]=-a; #Vector for the circulant matrix
        c[2]=1;
        c[N]=1;
        T=Circulant(c);
        T=convert(Matrix,T);
        k=nullspace(T);
        F=[T;k[:,1]';k[:,2]'];
        b1=dot(M2[1,2][2:N+1],k[:,1]);
        b2=dot(M2[1,2][2:N+1],k[:,2]);
        w=M2[1,1][1];
        c=(2-a)*M2[1,1][2:N+1];
        c=[c;b1;b2];
        x=F\c;
        b1=dot(M2[1,2][N+2:2*N+1],k[:,1]);
        b2=dot(M2[1,2][N+2:2*N+1],k[:,2]);
        c=(2-a)*M2[1,1][N+2:2*N+1];
        c=[c;b1;b2];
        y=F\c;
        M10new=[w;x;y];

        M10new=normalizemask(M10new);
        B[1:2*N+1,2]=M10new;

        M01=completeshift(M10new,N,1,"EVREGBORDER1");
        B[:,4]=M01;

        #Compute new M20
        M20new=(1/10)*(M2[1,6]-M2[1,1]+5*B[:,2]+10*M2[1,4]-5*M2[1,5]);

        M20new=normalizemask(M20new);
        B[:,3]=M20new;

        M02=completeshift(M20new,N,1,"EVREGBORDER1");
        B[:,7]=M02;
    else
        M20new=(1/10)*(M2[1,6]-M2[1,1]+5*M2[1,2]+10*M2[1,4]-5*M2[1,5]);

        M20new=normalizemask(M20new);
        B[:,3]=M20new;

        M02=completeshift(M20new,N,1,"EVREGBORDER1");
        B[:,7]=M02;
    end

    #M11

    if N%2==0   #Even case I modify here M20VER-->M11-->M30
        if circulant_system=="yes"
            v=zeros(N); #Circulant matrix vector
            k=zeros(1,N); #It will be the kernel vector -1 1 -1 ... 1
            for i in 1:N
                k[i]=(-1)^i;
            end
            v[1]=1; #Vector for the circulant matrix
            v[N]=1;
            A=convert(Matrix,Circulant(v)); #Circulant matrix for the M11 LS system
            Z=[A;k];

            #M20VER

            f=zeros(N);
            for i in 1:N
                f[i]=dot(B[2:N+1,3],A[i,:]);
            end
            ff=[f;0];
            xM20new=Z\ff;
            for i in 1:N
                f[i]=dot(B[N+2:2*N+1,3],A[i,:]);
            end
            ff=[f;0];
            yM20new=Z\ff;
            w=B[1,3];
            M20even_new=[w;xM20new;yM20new];

            M20even_new=normalizemask(M20even_new);
            B[1:2*N+1,3]=M20even_new;

            M02=completeshift(M20even_new,N,1,"EVREGBORDER1");
            B[:,7]=M02;

            #M11

            w=(1/10)*(a*M2[1,1][1]+(10-5*a)*B[1,2]+4*a*M20even_new[1]);
            c=(1/5)*(a*M2[1,1][2:N+1]+(10-5*a)*B[2:N+1,2]+4*a*B[2:N+1,3]);
            c=[c;dot(k,M2[2,2][2:N+1])];
            x=Z\c; #First neighborhood solution
            cc=(1/5)*(a*M2[1,1][N+2:2*N+1]+(10-5*a)*B[N+2:2*N+1,2]+4*a*B[N+2:2*N+1,3]);
            cc=[cc;dot(k,M2[2,2][N+2:2*N+1])];
            y=Z\cc; #Second neighborhood solution
            M11new=[w;x;y];

            M11new=normalizemask(M11new);
            B[1:2*N+1,5]=M11new;

            #M30

            M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);
            M30new=normalizemask(M30new);

            B[:,16]=M30new;

            M03new=completeshift(M30new,N,1,"EVREGBORDER1");

            B[:,30]=M03new;

        else

            #Compute M20 from M11 and then M30

            #M20

            M11_0=completeshift(M2[2,2],N,-1,"EVREGBORDER1");
            M20new=(1/(4*a))*(5*M2[2,2]-a*M2[1,1]+5*(a-2)*B[:,2]+5*M11_0);

            M20new=normalizemask(M20new);
            B[:,3]=M20new;

            M02=completeshift(M20new,N,1,"EVREGBORDER1");
            B[:,7]=M02;

            #M30

            M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);

            M30new=normalizemask(M30new);
            B[:,16]=M30new;

            M03new=completeshift(M30new,N,1,"EVREGBORDER1");

            M03new=normalizemask(M03new);
            B[:,30]=M03new;

        end

        elseif N%2!=0 #Odd case
            if circulant_system=="yes"

                #M11

                v=zeros(N);
                v[1]=1;
                v[N]=1;

                #Quadratic glueing data
                w=(1/10)*(a*M2[1,1][1]+(10-5*a)*B[1,2]+4*a*B[1,3]);
                A=convert(Matrix,Circulant(v));
                c=(1/5)*(a*M2[1,1][2:N+1]+(10-5*a)*B[2:N+1,2]+4*a*B[2:N+1,3]);
                x=A\c;
                c=(1/5)*(a*M2[1,1][N+2:2*N+1]+(10-5*a)*B[N+2:2*N+1,2]+4*a*B[N+2:2*N+1,3]);
                y=A\c;
                M11new=[w;x;y];

                M11new=normalizemask(M11new);
                B[1:2*N+1,5]=M11new;

            else

                #M20

                M11_0=completeshift(M2[2,2],N,-1,"EVREGBORDER1");
                M20new=(1/(4*a))*(5*M2[2,2]-a*M2[1,1]+5*(a-2)*B[:,2]+5*M11_0);

                M20new=normalizemask(M20new);
                B[:,3]=M20new;

                M02=completeshift(M20new,N,1,"EVREGBORDER1");
                B[:,7]=M02;

                #M30

                M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);

                M30new=normalizemask(M30new);
                B[:,16]=M30new;

                M03new=completeshift(M30new,N,1,"EVREGBORDER1");
                B[:,30]=M03new;

            end
        end

        #M21 non symmetric

        if symm=="no"

            M12_0=completeshift(M2[3,2],N,-1,"EVREGBORDER1");
            M21=(1/10)*(-a*M2[1,1]+5*a*B[:,2]+(20-10*a)*B[:,3]+6*a*B[:,16]-10*M12_0);

            M21=normalizemask(M21);
            B[:,6]=M21;

        #M21 symmetric

        else

            R=(1/10)*(-a*M2[1,1]+5*a*B[:,2]+(20-10*a)*B[:,3]+6*a*B[:,16]);
            A12_0=completeshift(M2[3,2],N,-1,"EVREGBORDER1");
            A21=M2[2,3];
            theta=(A21-A12_0);
            M21=1/2*(R+theta);

            M21=normalizemask(M21);
            B[:,6]=M21;

            M12_0=1/2*(R-theta);

            M12_0=normalizemask(M12_0);

            M12=completeshift(M12_0,N,1,"EVREGBORDER1");
            B[:,8]=M12;

        end

    #M31 non symmetric

    if symm=="no"

        M13_0=completeshift(M2[4,2],N,-1,"EVREGBORDER1");
        M31=(1/10)*(a*M2[1,1]-5*a*B[:,2]+10*a*B[:,3]+(20-10*a)*B[:,16]+4*a*M2[1,5]-10*M13_0);

        M31=normalizemask(M31);
        B[:,17]=M31;

    #M31 symmetric

    else

        R=(1/10)*(a*M2[1,1]-5*a*B[:,2]+10*a*B[:,3]+(20-10*a)*B[:,16]+4*a*M2[1,5]);
        A13_0=completeshift(M2[4,2],N,-1,"EVREGBORDER1");
        A31=M2[2,4];
        theta=(A31-A13_0);
        M31=1/2*(R+theta);

        M31=normalizemask(M31);
        B[:,17]=M31;

        M13_0=1/2*(R-theta);
        M13_0=normalizemask(M13_0);
        M13=completeshift(M13_0,N,1,"EVREGBORDER1");
        B[:,33]=M13;

    end

elseif T=="EVREGBORDER2"

    M2=degree_elevate_mask(N,"BORDER2",0); #7
    B=ordermatrix(M2,N,"EVREGBORDER2");
    a=2*cos(2*pi/N);

    if N!=3 #Compute the new M10

        c=zeros(N); #Right hand of the system
        c[1]=-a; #Vector for the circulant matrix
        c[2]=1;
        c[N]=1;
        T=Circulant(c);
        T=convert(Matrix,T);
        k=nullspace(T);
        F=[T;k[:,1]';k[:,2]'];
        b1=dot(M2[1,2][2:N+1],k[:,1]);
        b2=dot(M2[1,2][2:N+1],k[:,2]);
        w=M2[1,1][1];
        c=(2-a)*M2[1,1][2:N+1];
        c=[c;b1;b2];
        x=F\c;
        b1=dot(M2[1,2][N+2:2*N+1],k[:,1]);
        b2=dot(M2[1,2][N+2:2*N+1],k[:,2]);
        c=(2-a)*M2[1,1][N+2:2*N+1];
        c=[c;b1;b2];
        y=F\c;
        M10new=[w;x;y];

        M10new=normalizemask(M10new);
        B[1:2*N+1,2]=M10new;

        M01=completeshift(M10new,N,1,"EVREGBORDER2");
        B[:,4]=M01;

        #Compute new M20
        M20new=(1/10)*(M2[1,6]-M2[1,1]+5*B[:,2]+10*M2[1,4]-5*M2[1,5]);

        M20new=normalizemask(M20new);
        B[:,3]=M20new;

        M02=completeshift(M20new,N,1,"EVREGBORDER2");
        B[:,7]=M02;
    else
        M20new=(1/10)*(M2[1,6]-M2[1,1]+5*M2[1,2]+10*M2[1,4]-5*M2[1,5]);

        M20new=normalizemask(M20new);
        B[:,3]=M20new;

        M02=completeshift(M20new,N,1,"EVREGBORDER2");
        B[:,7]=M02;
    end

    #M11

    if N%2==0   #Even case I modify here M20VER-->M11-->M30
        if circulant_system=="yes"
            v=zeros(N); #Circulant matrix vector
            k=zeros(1,N); #It will be the kernel vector -1 1 -1 ... 1
            for i in 1:N
                k[i]=(-1)^i;
            end
            v[1]=1; #Vector for the circulant matrix
            v[N]=1;
            A=convert(Matrix,Circulant(v)); #Circulant matrix for the M11 LS system
            Z=[A;k];

            #M20VER

            f=zeros(N);
            for i in 1:N
                f[i]=dot(B[2:N+1,3],A[i,:]);
            end
            ff=[f;0];
            xM20new=Z\ff;
            for i in 1:N
                f[i]=dot(B[N+2:2*N+1,3],A[i,:]);
            end
            ff=[f;0];
            yM20new=Z\ff;
            w=B[1,3];
            M20even_new=[w;xM20new;yM20new];

            M20even_new=normalizemask(M20even_new);
            B[1:2*N+1,3]=M20even_new;

            M02=completeshift(M20even_new,N,1,"EVREGBORDER2");
            B[:,7]=M02;

            #M11

            w=(1/10)*(a*M2[1,1][1]+(10-5*a)*B[1,2]+4*a*M20even_new[1]);
            c=(1/5)*(a*M2[1,1][2:N+1]+(10-5*a)*B[2:N+1,2]+4*a*B[2:N+1,3]);
            c=[c;dot(k,M2[2,2][2:N+1])];
            x=Z\c; #First neighborhood solution
            cc=(1/5)*(a*M2[1,1][N+2:2*N+1]+(10-5*a)*B[N+2:2*N+1,2]+4*a*B[N+2:2*N+1,3]);
            cc=[cc;dot(k,M2[2,2][N+2:2*N+1])];
            y=Z\cc; #Second neighborhood solution
            M11new=[w;x;y];

            M11new=normalizemask(M11new);
            B[1:2*N+1,5]=M11new;

            #M30

            M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);
            M30new=normalizemask(M30new);

            B[:,16]=M30new;

            M03new=completeshift(M30new,N,1,"EVREGBORDER2");

            B[:,30]=M03new;

        else

            #Compute M20 from M11 and then M30

            #M20

            M11_0=completeshift(M2[2,2],N,-1,"EVREGBORDER2");
            M20new=(1/(4*a))*(5*M2[2,2]-a*M2[1,1]+5*(a-2)*B[:,2]+5*M11_0);

            M20new=normalizemask(M20new);
            B[:,3]=M20new;

            M02=completeshift(M20new,N,1,"EVREGBORDER2");
            B[:,7]=M02;

            #M30

            M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);

            M30new=normalizemask(M30new);
            B[:,16]=M30new;

            M03new=completeshift(M30new,N,1,"EVREGBORDER2");

            M03new=normalizemask(M03new);
            B[:,30]=M03new;

        end

        elseif N%2!=0 #Odd case
            if circulant_system=="yes"

                #M11

                v=zeros(N);
                v[1]=1;
                v[N]=1;

                #Quadratic glueing data
                w=(1/10)*(a*M2[1,1][1]+(10-5*a)*B[1,2]+4*a*B[1,3]);
                A=convert(Matrix,Circulant(v));
                c=(1/5)*(a*M2[1,1][2:N+1]+(10-5*a)*B[2:N+1,2]+4*a*B[2:N+1,3]);
                x=A\c;
                c=(1/5)*(a*M2[1,1][N+2:2*N+1]+(10-5*a)*B[N+2:2*N+1,2]+4*a*B[N+2:2*N+1,3]);
                y=A\c;
                M11new=[w;x;y];

                M11new=normalizemask(M11new);
                B[1:2*N+1,5]=M11new;

            else

                #M20

                M11_0=completeshift(M2[2,2],N,-1,"EVREGBORDER2");
                M20new=(1/(4*a))*(5*M2[2,2]-a*M2[1,1]+5*(a-2)*B[:,2]+5*M11_0);

                M20new=normalizemask(M20new);
                B[:,3]=M20new;

                M02=completeshift(M20new,N,1,"EVREGBORDER2");
                B[:,7]=M02;

                #M30

                M30new=(1/10)*(M2[1,1]-5*B[:,2]+10*B[:,3]+5*M2[1,5]-M2[1,6]);

                M30new=normalizemask(M30new);
                B[:,16]=M30new;

                M03new=completeshift(M30new,N,1,"EVREGBORDER2");
                B[:,30]=M03new;

            end
        end

        #M21 non symmetric

        if symm=="no"

            M12_0=completeshift(M2[3,2],N,-1,"EVREGBORDER2");
            M21=(1/10)*(-a*M2[1,1]+5*a*B[:,2]+(20-10*a)*B[:,3]+6*a*B[:,16]-10*M12_0);

            M21=normalizemask(M21);
            B[:,6]=M21;

        #M21 symmetric

        else

            R=(1/10)*(-a*M2[1,1]+5*a*B[:,2]+(20-10*a)*B[:,3]+6*a*B[:,16]);
            A12_0=completeshift(M2[3,2],N,-1,"EVREGBORDER2");
            A21=M2[2,3];
            theta=(A21-A12_0);
            M21=1/2*(R+theta);

            M21=normalizemask(M21);
            B[:,6]=M21;

            M12_0=1/2*(R-theta);

            M12_0=normalizemask(M12_0);

            M12=completeshift(M12_0,N,1,"EVREGBORDER2");
            B[:,8]=M12;

        end

    #M31 non symmetric

    if symm=="no"

        M13_0=completeshift(M2[4,2],N,-1,"EVREGBORDER2");
        M31=(1/10)*(a*M2[1,1]-5*a*B[:,2]+10*a*B[:,3]+(20-10*a)*B[:,16]+4*a*M2[1,5]-10*M13_0);

        M31=normalizemask(M31);
        B[:,17]=M31;

    #M31 symmetric

    else

        R=(1/10)*(a*M2[1,1]-5*a*B[:,2]+10*a*B[:,3]+(20-10*a)*B[:,16]+4*a*M2[1,5]);
        A13_0=completeshift(M2[4,2],N,-1,"EVREGBORDER2");
        A31=M2[2,4];
        theta=(A31-A13_0);
        M31=1/2*(R+theta);

        M31=normalizemask(M31);
        B[:,17]=M31;

        M13_0=1/2*(R-theta);
        M13_0=normalizemask(M13_0);
        M13=completeshift(M13_0,N,1,"EVREGBORDER2");
        B[:,33]=M13;

    end




else #T=="BORDEREV" for borders EVs we directly impose the G1 conditions

    a=2*cos(2*pi/N);
    k=convert(Int64,N/2);
    B1=degree_elevate_mask(N,"BORDEREVR",k); #3 #right border
    B2=degree_elevate_mask(N,"BORDEREVL",k); #4 #left border
    B3=degree_elevate_mask(N,"BORDERFACEEV",k); #5 #inner face border ev vertex
    B1=ordermatrix(B1,N,"BEV");
    B2=ordermatrix(B2,N,"BEV");
    B3=ordermatrix(B3,N,"BEV");
    Q=zeros(3*N+3,7*k-6);

    #M10 masks

    d=ones(k-1); #diagonal tridiag system
    l=1.0*ones(k-2); #lower diag
    u=l; #upper diag
    d=-a*d;
    A=convert(Matrix,Tridiagonal(l,d,u));
    M00=B1[:,1];
    MR=B1[:,2];
    ML=B2[:,4];
    ker=nullspace(A);
    A=[A;ker'];
    r=zeros(3*N+3,k);
    r[:,1]=(2-a)*M00-MR;
    r[:,2:end-2].=(2-a)*M00;
    r[:,end-1]=(2-a)*M00-ML;
    M10_ACC5=B1[:,4];
    r[:,end]=r[:,end]+M10_ACC5*ker[1];
    M10_ACC5=B3[:,4];

    for i in 2:k-2
        r[:,end]=r[:,end]+M10_ACC5*ker[i];
        M10_ACC5=completeshift(M10_ACC5,k,1,"BORDER");
    end
    M10_ACC5=B2[:,2];
    r[:,end]=r[:,end]+M10_ACC5*ker[end];

    for i in 1:3*N+3
        sol=A\r[i,:];
        Q[i,1:k-1]=sol;
    end

    for i in 1:k-1
        Q[:,i]=normalizemask(Q[:,i]);
    end

    #M20

    M30=B3[:,16];
    M40=B3[:,13];
    M50=B3[:,10];

    for i in k:2*k-2
        Q[:,i]=(1/10)*(-M00+5*Q[:,i-k+1]+10*M30-5*M40+M50);
        M30=completeshift(M30,k,1,"BORDER");
        M40=completeshift(M40,k,1,"BORDER");
        M50=completeshift(M50,k,1,"BORDER");
    end

    for i in 3:4
        Q[:,i]=normalizemask(Q[:,i]);
    end

    #M11

    C=zeros(k-1,k);
    v=zeros(1,k);
    v[1]=1;
    v[2]=1;
    C[1,:]=v;
    v=vectorshift(v,1);

    for i in 2:k-1
        C[i,:]=v;
        v=vectorshift(v,1);
    end

    r=zeros(3*N+3,k-1);
    s=zeros(3*N+3,k-1);
    H=zeros(3*N+3,k);
    H[:,1]=B1[:,5];
    H[:,end]=B2[:,5];
    M11_ACC5=B3[:,5];

    for i in 2:k-1
        H[:,i]=M11_ACC5;
        M11_ACC5=completeshift(M11_ACC5,k,1,"BORDER");
    end

    for i in 1:k-1
        r[:,i]=(1/5)*(a*M00+5*(2-a)*Q[:,i]+4*a*Q[:,i+k-1]);
    end

    for i in 1:3*N+3
        s[i,:]=-C*H[i,:];
    end

    r=r+s;

    for i in 1:3*N+3
        sol=(C*transpose(C))\r[i,:]; #compute the lambdas of the projection
        Q[i,2*k-1:3*k-2]=H[i,:]+transpose(C)*sol;
    end

    for i in 2*k-1:3*k-2
        Q[:,i]=normalizemask(Q[:,i]);
    end

    #M21-M12

    M30=B3[:,16];
    M12=B1[:,8];
    M21=B3[:,6];
    R=(1/10)*(-a*M00+5*a*Q[:,1]-10*(a-2)*Q[:,k]+6*a*M30);
    theta=M21-M12;
    Q[:,3*k-1]=(1/2)*(R-theta);
    Q[:,4*k-2]=(1/2)*(R+theta);
    M12=B3[:,8];

    for i in 2:k-2
        M30=completeshift(M30,k,1,"BORDER");
        R=(1/10)*(-a*M00+5*a*Q[:,i]-10*(a-2)*Q[:,k+i-1]+6*a*M30);
        M21=completeshift(M21,k,1,"BORDER");
        theta=M21-M12;
        Q[:,3*k+i-2]=(1/2)*(R-theta);
        Q[:,4*k+i-3]=(1/2)*(R+theta);
        M12=completeshift(M12,k,1,"BORDER");
    end

    M21=B2[:,6];
    M30=B2[:,16];
    R=(1/10)*(-a*M00+5*a*Q[:,k-1]-10*(a-2)*Q[:,2*k-2]+6*a*M30);
    theta=M21-M12;
    Q[:,4*k-3]=(1/2)*(R-theta);
    Q[:,5*k-4]=(1/2)*(R+theta);

    for i in 3*k-1:5*k-4
        Q[:,i]=normalizemask(Q[:,i]);
    end

    #M31-M13

    M40=B3[:,13];
    M30=B3[:,16];
    M13=B1[:,33];
    M31=B3[:,17];
    R2=(1/10)*(a*M00-5*a*Q[:,1]+10*a*Q[:,k]-10*(a-2)*M30+4*a*M40);
    theta2=M31-M13;
    Q[:,5*k-3]=(1/2)*(R2-theta2);
    Q[:,6*k-4]=(1/2)*(R2+theta2);
    M13=B3[:,33];

    for i in 2:k-2
        M30=completeshift(M30,k,1,"BORDER");
        M40=completeshift(M40,k,1,"BORDER");
        R2=(1/10)*(a*M00-5*a*Q[:,i]+10*a*Q[:,k+i-1]-10*(a-2)*M30+4*a*M40);
        M31=completeshift(M31,k,1,"BORDER");
        theta2=M31-M13;
        Q[:,5*k+i-4]=(1/2)*(R2-theta2);
        Q[:,6*k+i-5]=(1/2)*(R2+theta2);
        M13=completeshift(M13,k,1,"BORDER");
    end

    M31=B2[:,17];
    M30=B2[:,16];
    M40=B2[:,13];
    R2=(1/10)*(a*M00-5*a*Q[:,k-1]+10*a*Q[:,2*k-2]-10*(a-2)*M30+4*a*M40);
    theta2=M31-M13;
    Q[:,6*k-5]=(1/2)*(R2-theta2);
    Q[:,7*k-6]=(1/2)*(R2+theta2);

    for i in 5*k-3:7*k-6
        Q[:,i]=normalizemask(Q[:,i]);
    end

B=Q;

end

return B;

end
