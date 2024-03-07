#Basis degree elevation/knots insertion

currentdir =  @__DIR__

include(joinpath(currentdir,"G1Splines.jl"));

meshname="triangle_planar";

if (size(ARGS,1)>0)
    m = offread(ARGS[1]);
else
    m = offread(joinpath(currentdir,string("../data/",meshname,".off")));
end

hm = hmesh(m);

N=3;

#M=g1basis_EV(N);

knt=knots(5,1,1)
mm, d = dim_deg(knt)
kntelev=knots(5,2,3)
melev, delev = dim_deg(kntelev)

resold=acc_d(hm,d); #Collection of spline patches #res=[sup,gs]
gsurfold=resold;
res=acc_d(hm,melev-1); #Collection of spline patches #res=[sup,gs]
gsurf=res;

basis=assemble_g1basis_spline(hm,kntelev);

if true
    for i in 1:size(basis)[2]
        println(size(basis)[2])
        row=findnz(basis[:,i])[1];
        val=findnz(basis[:,i])[2];
        local gsbas = GSpline(0, kntelev,hm )
        CP=copy(gsurf.ctrpoints);
        for j in 1:length(row)
            CP[3,row[j]]=val[j];
        end
        gsbas.ctrpoints=CP;
        @axlview gsbas







    end
end

#BB=bezier_to_spline(M,kntelev,N);



#=for i in 1:size(M,2)

    b=M[:,i];
    data=findnz(b);
    row=data[1];
    val=data[2];
    B=zeros(m,m,N);

    for j=1:length(row)
        p=div(row[j],m^2)+1;
        println(p)
        h=div(row[j]%m^2,m)+1;
        w=row[j]%m;
        B[w,h,p]=val[j];
    end

    #=oldctps=zeros(1,0);

    for i in 1:N
        for j in 1:m
            oldctps=hcat(oldctps,B[j,:,i]');
        end
    end

    local gsbas = GSpline(0, knt,hm )
        CP=copy(gsurfold.ctrpoints);
        CP[3,:]=oldctps;
        gsbas.ctrpoints=CP;
        @axlview gsbas=#

    Belev=zeros(m+1,m+1,N)

    for i in 1:N
        Belev[:,:,i]=degree_elevate(B[:,:,i]);
    end

    newcpts=zeros(0);

    for i in 1:N
        for j in 1:m+1
            newcpts=vcat(newcpts,Belev[j,:,i]);
        end
    end

    newcpts=sparse(newcpts)

          row=findnz(newcpts)[1];
          val=findnz(newcpts)[2];
          local gsbas = GSpline(0, kntelev,hm )
          CP=copy(gsurf.ctrpoints);
          for j in 1:length(row)
              CP[3,row[j]]=val[j];
          end
          gsbas.ctrpoints=CP;
          @axlview gsbas

    #=local gsbas = GSpline(0, kntelev,hm )
        CP=copy(gsurf.ctrpoints);
        CP[3,:]=newcpts;
        gsbas.ctrpoints=CP;
        #@axlview gsbas=#



end=#
