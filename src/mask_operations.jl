function vectorshift(v::Array, s::Int) #This function return a vector shifted forward if s==1, backward if s==-1

    L=length(v);
    h=zeros(L);
    if s==1 #forward
        h[1]=v[L];
        for i in 2:L
            h[i]=v[i-1];
        end
    elseif s== -1 #backward
        h[L]=v[1];
        for i in 1:L-1
            h[i]=v[i+1];
        end
    end
    return h;
end


function completeshift(v::Array,N::Int64,p::Int64,S::String="BORDER") #This function operates a complete shift of a subdivision mask

  if S=="BORDER"
    u=zeros(6*N+3,1);
    u[1]=v[1];
    u[2:N+2]=vectorshift(v[2:N+2],p);
    u[N+3:2*N+2]=vectorshift(v[N+3:2*N+2],p);
    u[2*N+3:3*N+3]=vectorshift(v[2*N+3:3*N+3],p);
    u[3*N+4:4*N+3]=vectorshift(v[3*N+4:4*N+3],p);
    u[4*N+4:5*N+3]=vectorshift(v[4*N+4:5*N+3],p);
    u[5*N+4:6*N+3]=vectorshift(v[5*N+4:6*N+3],p);

  elseif S=="EVREGBORDER"
    u=zeros(2*N+1,1);
    u[1]=v[1];
    u[2:N+1]=vectorshift(v[2:N+1],p);
    u[N+2:2*N+1]=vectorshift(v[N+2:2*N+1],p);

  else
  u=zeros(6*N+1,1);
  if length(v)<length(u)
      u[1]=v[1];
      u[2:N+1]=vectorshift(v[2:N+1],p);
      u[N+2:2*N+1]=vectorshift(v[N+2:2*N+1],p);
  else
    u[1]=v[1];
    u[2:N+1]=vectorshift(v[2:N+1],p);
    u[N+2:2*N+1]=vectorshift(v[N+2:2*N+1],p);
    u[2*N+2:3*N+1]=vectorshift(v[2*N+2:3*N+1],p);
    u[3*N+2:4*N+1]=vectorshift(v[3*N+2:4*N+1],p);
    u[4*N+2:5*N+1]=vectorshift(v[4*N+2:5*N+1],p);
    u[5*N+2:6*N+1]=vectorshift(v[5*N+2:6*N+1],p);
  end
  end

  return u

  end

  function normalizemask(v::Array) #This function normalize the weights of the subdivision masks
  norm=sum(v);
  if abs(norm)<=1.e-10
      norm=1;
  end
  v=v/norm;

  return v

  end
