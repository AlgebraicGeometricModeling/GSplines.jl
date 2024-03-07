function acc3bigmatrix(N::Int, S::String)
    #N valence of the vertex
    if S=="INNEREV"
        d=3; #cubic surface
        row=6*N+1; #col is the numbers of column of M
        row=convert(Int64,row);
        col=(d+1)^2; #col is the numbers of column of M
        col=convert(Int64,col);
        M=zeros(row,col);
        rho=1/(N*(N+5));

        M[1,1]=rho*N^2;
        M[2:N+1,1].=4*rho;
        M[N+2:2*N+1,1].=rho;

        M[1,2]=rho*N^2;
        M[2,2]=2*rho*N;
        M[3,2]=rho*N;
        M[N+1,2]=rho*N;
        M[N+2,2]=(rho*N)/2;
        M[2*N+1,2]=(rho*N)/2;

        M[1,3]=rho*N^2;
        M[2,3]=rho*N;
        M[3,3]=2*rho*N;
        M[4,3]=rho*N;
        M[N+2,3]=(rho*N)/2;
        M[N+3,3]=(rho*N)/2;

        M[1,4]=rho*N^2;
        M[2,4]=2*rho*N;
        M[3,4]=2*rho*N;
        M[N+2,4]=rho*N;

        M[1,5]=1/9;
        M[2,5]=4/9;
        M[3,5]=1/36;
        M[N+1,5]=1/36;
        M[N+2,5]=1/9;
        M[2*N+1,5]=1/9;
        M[2*N+2,5]=1/9;
        M[3*N+2,5]=1/36;
        M[end,5]=1/36;

        M[1,6]=1/9;
        M[2,6]=4/9;
        M[3,6]=1/18;
        M[N+2,6]=2/9;
        M[2*N+2,6]=1/9;
        M[3*N+2,6]=1/18;

        M[1,7]=2/9;
        M[2,7]=4/9;
        M[3,7]=1/18;
        M[N+1,7]=1/18;
        M[N+2,7]=1/9;
        M[2*N+1,7]=1/9;

        M[1,8]=2/9;
        M[2,8]=4/9;
        M[3,8]=1/9;
        M[N+2,8]=2/9;

        M[1,9]=1/36;
        M[2,9]=1/9;
        M[3,9]=1/9;
        M[N+2,9]=4/9;
        M[2*N+2,9]=1/36;
        M[2*N+3,9]=1/36;
        M[3*N+2,9]=1/9;
        M[4*N+2,9]=1/36;
        M[5*N+2,9]=1/9;

        M[1,10]=1/18;
        M[2,10]=1/9;
        M[3,10]=2/9;
        M[N+2,10]=4/9;
        M[2*N+3,10]=1/18;
        M[5*N+2,10]=1/9;

        M[1,11]=1/18;
        M[2,11]=2/9;
        M[3,11]=1/9;
        M[N+2,11]=4/9;
        M[2*N+2,11]=1/18;
        M[3*N+2,11]=1/9;

        M[1,12]=1/9;
        M[2,12]=2/9;
        M[3,12]=2/9;
        M[N+2,12]=4/9;

        M[1,13]=1/9;
        M[2,13]=1/36;
        M[3,13]=4/9;
        M[4,13]=1/36;
        M[N+2,13]=1/9;
        M[N+3,13]=1/9;
        M[2*N+3,13]=1/9;
        M[3*N+3,13]=1/36;
        M[5*N+2,13]=1/36;

        M[1,14]=2/9;
        M[2,14]=1/18;
        M[3,14]=4/9;
        M[4,14]=1/18;
        M[N+2,14]=1/9;
        M[N+3,14]=1/9;

        M[1,15]=1/9;
        M[2,15]=1/18;
        M[3,15]=4/9;
        M[N+2,15]=2/9;
        M[2*N+3,15]=1/9;
        M[5*N+2,15]=1/18;

        M[1,16]=2/9;
        M[2,16]=1/9;
        M[3,16]=4/9;
        M[N+2,16]=2/9;

    elseif S=="BORDER1"

        d=3; #cubic surface
        row=2*N+4; #row is the numbers of rows of M
        row=convert(Int64,row);
        col=(d+1)^2; #col is the numbers of columns of M
        col=convert(Int64,col);
        M=zeros(row,col);
        k=2;#convert(Int64,N/2);
        rho=1/(N*(N+5));
        M[1,1]=rho*N^2;
        M[2:N+1,1].=4*rho;
        M[N+2:2*N+1,1].=rho;

        M[1,2]=rho*N^2;
        M[2,2]=2*rho*N;
        M[3,2]=rho*N;
        M[N+1,2]=rho*N;
        M[N+2,2]=(rho*N)/2;
        M[2*N+1,2]=(rho*N)/2;

        M[1,3]=rho*N^2;
        M[2,3]=rho*N;
        M[3,3]=2*rho*N;
        M[4,3]=rho*N;
        M[N+2,3]=(rho*N)/2;
        M[N+3,3]=(rho*N)/2;

        M[1,4]=rho*N^2;
        M[2,4]=2*rho*N;
        M[3,4]=2*rho*N;
        M[N+2,4]=rho*N;

        M[1,5]=1/9;
        M[2,5]=4/9;
        M[3,5]=1/36;
        M[N+1,5]=1/36;
        M[N+2,5]=1/9;
        M[2*N+1,5]=1/9;
        M[2*N+2,5]=1/9;
        M[2*N+3,5]=1/36;
        M[2*N+4,5]=1/36;

        M[1,6]=1/9;
        M[2,6]=4/9;
        M[3,6]=1/18;
        M[N+2,6]=2/9;
        M[2*N+2,6]=1/9;
        M[2*N+3,6]=1/18;

        M[1,7]=2/9;
        M[2,7]=4/9;
        M[3,7]=1/18;
        M[N+1,7]=1/18;
        M[N+2,7]=1/9;
        M[2*N+1,7]=1/9;

        M[1,8]=2/9;
        M[2,8]=4/9;
        M[3,8]=1/9;
        M[N+2,8]=2/9;

        M[3,9]=1/6;
        M[6,9]=2/3;
        M[2*N+3,9]=1/6;

        M[3,10]=1/3;
        M[6,10]=2/3;

        M[1,11]=1/18;
        M[2,11]=2/9;
        M[3,11]=1/9;
        M[N+2,11]=4/9;
        M[2*N+2,11]=1/18;
        M[2*N+3,11]=1/9;

        M[1,12]=1/(2*k+5);
        M[2,12]=2/(2*k+5);
        M[3,12]=2/(2*k+5);
        M[N+2,12]=2*k/(2*k+5);

        M[3,13]=2/3;
        M[N+2,13]=1/6;
        M[N+3,13]=1/6;


        M[1,14]=2/9;
        M[2,14]=1/18;
        M[3,14]=4/9;
        M[4,14]=1/18;
        M[N+2,14]=1/9;
        M[N+3,14]=1/9;


        M[3,15]=2/3;
        M[N+2,15]=1/3;


        M[1,16]=2/(2*k+5);
        M[2,16]=1/(2*k+5);
        M[3,16]=2*k/(2*k+5);
        M[N+2,16]=2/(2*k+5);

    elseif S=="BORDER2"

        d=3; #cubic surface
        row=2*N+4; #row is the numbers of rows of M
        row=convert(Int64,row);
        col=(d+1)^2; #col is the numbers of columns of M
        col=convert(Int64,col);
        M=zeros(row,col);
        k=2;#convert(Int64,N/2);
        rho=1/(N*(N+5));
        M[1,1]=rho*N^2;
        M[2:N+1,1].=4*rho;
        M[N+2:2*N+1,1].=rho;

        M[1,2]=rho*N^2;
        M[2,2]=2*rho*N;
        M[3,2]=rho*N;
        M[N+1,2]=rho*N;
        M[N+2,2]=(rho*N)/2;
        M[2*N+1,2]=(rho*N)/2;

        M[1,3]=rho*N^2;
        M[2,3]=rho*N;
        M[3,3]=2*rho*N;
        M[4,3]=rho*N;
        M[N+2,3]=(rho*N)/2;
        M[N+3,3]=(rho*N)/2;

        M[1,4]=rho*N^2;
        M[2,4]=2*rho*N;
        M[3,4]=2*rho*N;
        M[N+2,4]=rho*N;

        M[N+2,5]=1/6;
        M[2,5]=2/3;
        M[2*N+1,5]=1/6;

        M[N+2,6]=1/3;
        M[2,6]=2/3;

        M[1,7]=2/9;
        M[2,7]=4/9;
        M[3,7]=1/18;
        M[N+1,7]=1/18;
        M[N+2,7]=1/9;
        M[2*N+1,7]=1/9;

        M[3,8]=1/(2*k+5);
        M[N+2,8]=2/(2*k+5);
        M[1,8]=2/(2*k+5);
        M[2,8]=2*k/(2*k+5);


        M[N+2,9]=2/3;
        M[2,9]=1/6;
        M[2*N+3,9]=1/6;


        M[1,10]=1/18;
        M[2,10]=1/9;
        M[3,10]=2/9;
        M[N+2,10]=4/9;
        M[2*N+2,10]=1/18;
        M[2*N+3,10]=1/9;


        M[N+2,11]=2/3;
        M[2,11]=1/3;


        M[2,12]=2/(2*k+5);
        M[1,12]=1/(2*k+5);
        M[N+2,12]=2*k/(2*k+5);
        M[3,12]=2/(2*k+5);


        M[1,13]=1/9;
        M[3,13]=4/9;
        M[2,13]=1/36;
        M[4,13]=1/36;
        M[N+2,13]=1/9;
        M[N+3,13]=1/9;
        M[2*N+2,13]=1/9;
        M[2*N+3,13]=1/36;
        M[2*N+4,13]=1/36;

        M[N+2,14]=1/9;
        M[3,14]=4/9;
        M[2,14]=1/18;
        M[1,14]=2/9;
        M[N+3,14]=1/9;
        M[4,14]=1/18;

        M[N+2,15]=2/9;
        M[3,15]=4/9;
        M[2,15]=1/18;
        M[2*N+3,15]=1/18;
        M[1,15]=1/9;
        M[2*N+2,15]=1/9;

        M[1,16]=2/9;
        M[3,16]=4/9;
        M[2,16]=1/9;
        M[N+2,16]=2/9;


    elseif S=="CORNER"

        d=3; #cubic surface
        row=2*N+1; #row is the numbers of rows of M
        row=convert(Int64,row);
        k=2;
        col=(d+1)^2; #col is the numbers of columns of M
        col=convert(Int64,col);
        M=zeros(row,col);
        rho=1/(N*(N+5));

        M[1,1]=rho*N^2;
        M[2:N+1,1].=4*rho;
        M[N+2:2*N+1,1].=rho;

        M[1,2]=rho*N^2;
        M[2,2]=2*rho*N;
        M[3,2]=rho*N;
        M[N+1,2]=rho*N;
        M[N+2,2]=(rho*N)/2;
        M[2*N+1,2]=(rho*N)/2;

        M[1,3]=rho*N^2;
        M[2,3]=rho*N;
        M[3,3]=2*rho*N;
        M[4,3]=rho*N;
        M[N+2,3]=(rho*N)/2;
        M[N+3,3]=(rho*N)/2;

        M[1,4]=rho*N^2;
        M[2,4]=2*rho*N;
        M[3,4]=2*rho*N;
        M[N+2,4]=rho*N;

        M[2,5]=2/3;
        M[N+2,5]=1/6;
        M[2*N+1,5]=1/6;

        M[2,6]=2/3;
        M[N+2,6]=1/3;

        M[1,7]=2/9;
        M[2,7]=4/9;
        M[3,7]=1/18;
        M[N+1,7]=1/18;
        M[N+2,7]=1/9;
        M[2*N+1,7]=1/9;

        M[1,8]=2/(2*k+5);
        M[2,8]=2*k/(2*k+5);
        M[3,8]=1/(2*k+5);
        M[N+2,8]=2/(2*k+5);


        M[N+2,9]=1;


        M[3,10]=1/3;
        M[N+2,10]=2/3;


        M[2,11]=1/3;
        M[N+2,11]=2/3;


        M[1,12]=1/9;
        M[2,12]=2/9;
        M[3,12]=2/9;
        M[N+2,12]=4/9;

        M[3,13]=2/3;
        M[N+2,13]=1/6;
        M[N+3,13]=1/6;


        M[1,14]=2/9;
        M[2,14]=1/18;
        M[3,14]=4/9;
        M[4,14]=1/18;
        M[N+2,14]=1/9;
        M[N+3,14]=1/9;


        M[3,15]=2/3;
        M[N+2,15]=1/3;


        M[1,16]=2/(2*k+5);
        M[2,16]=1/(2*k+5);
        M[3,16]=2*k/(2*k+5);
        M[N+2,16]=2/(2*k+5);

    elseif S=="BORDEREVL"

        d=3; #cubic surface
        k=convert(Int64,N/2);
        row=3*N+3; #col is the numbers of column of M
        row=convert(Int64,row);
        col=(d+1)^2; #col is the numbers of column of M
        col=convert(Int64,col);
        M=zeros(row,col);
        rho=1/(N*(N+5));

        M[1,1]=2/3;
        M[2,1]=1/6;
        M[k+2,1]=1/6;

        M[1,2]=rho*N^2;
        M[k,2]=rho*N;
        M[k+1,2]=2*rho*N;
        M[k+2,2]=rho*N;
        M[2*k+1,2]=(rho*N)/2;
        M[2*k+2,2]=(rho*N)/2;

        M[1,3]=2/3;
        M[k+2,3]=1/3;

        M[1,4]=N/(N+5);
        M[k+1,4]=2/(N+5);
        M[k+2,4]=2/(N+5);
        M[2*k+2,4]=1/(N+5);

        M[1,5]=1/9;
        M[k,5]=1/36;
        M[k+1,5]=4/9;
        M[k+2,5]=1/36;
        M[2*k+1,5]=1/9;
        M[2*k+2,5]=1/9;
        M[3*k+2,5]=1/9;
        M[4*k+3,5]=1/36;
        M[6*k+2,5]=1/36;

        M[1,6]=1/9;
        M[k+1,6]=4/9;
        M[k+2,6]=1/18;
        M[2*k+2,6]=2/9;
        M[3*k+2,6]=1/9;
        M[4*k+3,6]=1/18;

        M[1,7]=2/9;
        M[k,7]=1/18;
        M[k+1,7]=4/9;
        M[k+2,7]=1/18;
        M[2*k+1,7]=1/9;
        M[2*k+2,7]=1/9;

        M[1,8]=2/9;
        M[k+1,8]=4/9;
        M[k+2,8]=1/9;
        M[2*k+2,8]=2/9;

        M[1,9]=1/36;
        M[k+1,9]=1/9;
        M[k+2,9]=1/9;
        M[2*k+2,9]=4/9;
        M[3*k+2,9]=1/36;
        M[3*k+3,9]=1/36;
        M[4*k+3,9]=1/9;
        M[5*k+3,9]=1/36;
        M[6*k+3,9]=1/9;

        M[1,10]=1/18;
        M[k+1,10]=1/9;
        M[k+2,10]=2/9;
        M[2*k+2,10]=4/9;
        M[3*k+3,10]=1/18;
        M[6*k+3,10]=1/9;

        M[1,11]=1/18;
        M[k+1,11]=2/9;
        M[k+2,11]=1/9;
        M[2*k+2,11]=4/9;
        M[3*k+2,11]=1/18;
        M[4*k+3,11]=1/9;

        M[1,12]=1/9;
        M[k+1,12]=2/9;
        M[k+2,12]=2/9;
        M[2*k+2,12]=4/9;

        M[1,13]=1/6;
        #M[2,13]=1/36;
        M[k+2,13]=2/3;
        #M[4,13]=1/36;
        #M[N+2,13]=1/9;
        #M[N+3,13]=1/9;
        M[3*k+3,13]=1/6;
        #M[3*N+3,13]=1/36;
        #M[5*N+2,13]=1/36;

        M[1,14]=1/3;
        M[k+2,14]=2/3;

        M[1,15]=1/9;
        M[k+1,15]=1/18;
        M[k+2,15]=4/9;
        M[2*k+2,15]=2/9;
        M[3*k+3,15]=1/9;
        M[6*k+3,15]=1/18;

        M[1,16]=2/9;
        M[k+1,16]=1/9;
        M[k+2,16]=4/9;
        M[2*k+2,16]=2/9;

    elseif S=="BORDEREVR"

        d=3; #cubic surface
        k=convert(Int64,N/2);
        row=3*N+3; #col is the numbers of column of M
        row=convert(Int64,row);
        col=(d+1)^2; #col is the numbers of column of M
        col=convert(Int64,col);
        M=zeros(row,col);
        rho=1/(N*(N+5));

        M[1,1]=2/3;
        M[2,1]=1/6;
        M[k+2,1]=1/6;

        M[1,2]=2/3;
        M[2,2]=1/3;

        M[1,3]=rho*N^2;
        M[2,3]=rho*N;
        M[3,3]=2*rho*N;
        M[4,3]=rho*N;
        M[k+3,3]=(rho*N)/2;
        M[k+4,3]=(rho*N)/2;

        M[1,4]=N/(N+5);
        M[2,4]=2/(N+5);
        M[3,4]=2/(N+5);
        M[k+3,4]=1/(N+5);

        M[1,5]=1/6;
        M[2,5]=2/3;
        #M[3,5]=1/36;
        #M[N+1,5]=1/36;
        #M[N+2,5]=1/9;
        #M[2*N+1,5]=1/9;
        M[2*k+3,5]=1/6;
        #M[3*N+2,5]=1/36;
        #M[end,5]=1/36;

        M[1,6]=1/9;
        M[2,6]=4/9;
        M[3,6]=1/18;
        M[k+3,6]=2/9;
        M[2*k+3,6]=1/9;
        M[3*k+4,6]=1/18;

        M[1,7]=1/3;
        M[2,7]=2/3;

        M[1,8]=2/9;
        M[2,8]=4/9;
        M[3,8]=1/9;
        M[k+3,8]=2/9;

        M[1,9]=1/36;
        M[2,9]=1/9;
        M[3,9]=1/9;
        M[k+3,9]=4/9;
        M[2*k+3,9]=1/36;
        M[2*k+4,9]=1/36;
        M[3*k+4,9]=1/9;
        M[4*k+4,9]=1/36;
        M[5*k+4,9]=1/9;

        M[1,10]=1/18;
        M[2,10]=1/9;
        M[3,10]=2/9;
        M[k+3,10]=4/9;
        M[2*k+4,10]=1/18;
        M[5*k+4,10]=1/9;

        M[1,11]=1/18;
        M[2,11]=2/9;
        M[3,11]=1/9;
        M[k+3,11]=4/9;
        M[2*k+3,11]=1/18;
        M[3*k+4,11]=1/9;

        M[1,12]=1/9;
        M[2,12]=2/9;
        M[3,12]=2/9;
        M[k+3,12]=4/9;

        M[1,13]=1/9;
        M[2,13]=1/36;
        M[3,13]=4/9;
        M[4,13]=1/36;
        M[k+3,13]=1/9;
        M[k+4,13]=1/9;
        M[2*k+4,13]=1/9;
        M[3*k+5,13]=1/36;
        M[5*k+4,13]=1/36;

        M[1,14]=2/9;
        M[2,14]=1/18;
        M[3,14]=4/9;
        M[4,14]=1/18;
        M[k+3,14]=1/9;
        M[k+4,14]=1/9;

        M[1,15]=1/9;
        M[2,15]=1/18;
        M[3,15]=4/9;
        M[k+3,15]=2/9;
        M[2*k+4,15]=1/9;
        M[5*k+4,15]=1/18;

        M[1,16]=2/9;
        M[2,16]=1/9;
        M[3,16]=4/9;
        M[k+3,16]=2/9;

    elseif S=="BORDERFACEEV"

        d=3; #cubic surface
        k=convert(Int64,N/2);
        row=3*N+3; #col is the numbers of column of M
        row=convert(Int64,row);
        col=(d+1)^2; #col is the numbers of column of M
        col=convert(Int64,col);
        M=zeros(row,col);
        rho=1/(N*(N+5));

        M[1,1]=2/3;
        M[2,1]=1/6;
        M[k+2,1]=1/6;

        M[1,2]=rho*N^2;
        M[2,2]=rho*N;
        M[3,2]=2*rho*N;
        M[4,2]=rho*N;
        M[k+3,2]=(rho*N)/2;
        M[k+4,2]=(rho*N)/2;

        M[1,3]=rho*N^2;
        M[3,3]=rho*N;
        M[4,3]=2*rho*N;
        M[5,3]=rho*N;
        M[k+4,3]=(rho*N)/2;
        M[k+5,3]=(rho*N)/2;

        M[1,4]=N/(N+5);
        M[3,4]=2/(N+5);
        M[4,4]=2/(N+5);
        M[k+4,4]=1/(N+5);

        M[1,5]=1/9;
        M[2,5]=1/36;
        M[3,5]=4/9;
        M[4,5]=1/36;
        M[k+3,5]=1/9;
        M[k+4,5]=1/9;
        M[2*k+4,5]=1/9;
        M[3*k+5,5]=1/36;
        M[5*k+4,5]=1/36;

        M[1,6]=1/9;
        M[3,6]=4/9;
        M[4,6]=1/18;
        M[k+4,6]=2/9;
        M[2*k+4,6]=1/9;
        M[3*k+5,6]=1/18;

        M[1,7]=2/9;
        M[2,7]=1/18;
        M[3,7]=4/9;
        M[4,7]=1/18;
        M[k+3,7]=1/9;
        M[k+4,7]=1/9;

        M[1,8]=2/9;
        M[3,8]=4/9;
        M[4,8]=1/9;
        M[k+4,8]=2/9;

        M[1,9]=1/36;
        M[3,9]=1/9;
        M[4,9]=1/9;
        M[k+4,9]=4/9;
        M[2*k+4,9]=1/36;
        M[2*k+5,9]=1/36;
        M[3*k+5,9]=1/9;
        M[4*k+5,9]=1/36;
        M[5*k+5,9]=1/9;

        M[1,10]=1/18;
        M[3,10]=1/9;
        M[4,10]=2/9;
        M[k+4,10]=4/9;
        M[2*k+5,10]=1/18;
        M[5*k+5,10]=1/9;

        M[1,11]=1/18;
        M[3,11]=2/9;
        M[4,11]=1/9;
        M[k+4,11]=4/9;
        M[2*k+4,11]=1/18;
        M[3*k+5,11]=1/9;

        M[1,12]=1/9;
        M[3,12]=2/9;
        M[4,12]=2/9;
        M[k+4,12]=4/9;

        M[1,13]=1/9;
        M[3,13]=1/36;
        M[4,13]=4/9;
        M[5,13]=1/36;
        M[k+4,13]=1/9;
        M[k+5,13]=1/9;
        M[2*k+5,13]=1/9;
        M[3*k+6,13]=1/36;
        M[5*k+5,13]=1/36;

        M[1,14]=2/9;
        M[3,14]=1/18;
        M[4,14]=4/9;
        M[5,14]=1/18;
        M[k+4,14]=1/9;
        M[k+5,14]=1/9;

        M[1,15]=1/9;
        M[3,15]=1/18;
        M[4,15]=4/9;
        M[k+4,15]=2/9;
        M[2*k+5,15]=1/9;
        M[5*k+5,15]=1/18;

        M[1,16]=2/9;
        M[3,16]=1/9;
        M[4,16]=4/9;
        M[k+4,16]=2/9;

    end

    return M;
end
