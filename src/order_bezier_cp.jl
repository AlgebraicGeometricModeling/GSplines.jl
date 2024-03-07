function order_bezier_cp(gs::GSpline,C::Array,P::Array,e::Int64,deg::Int64)
    #m=gs.mesh;
if deg==5
    #I start ordering from e

    C[:,idx(gs,e,1,1)]=P[:,1] #b00
    C[:,idx(gs,e,2,1)]=P[:,2] #b10
    C[:,idx(gs,e,3,1)]=P[:,3] #b20
    C[:,idx(gs,e,4,1)]=P[:,16] #b30
    C[:,idx(gs,e,5,1)]=P[:,13] #b40
    C[:,idx(gs,e,6,1)]=P[:,10] #b50

    C[:,idx(gs,e,1,2)]=P[:,4] #b01
    C[:,idx(gs,e,2,2)]=P[:,5] #b11
    C[:,idx(gs,e,3,2)]=P[:,6] #b21
    C[:,idx(gs,e,4,2)]=P[:,17] #b31
    C[:,idx(gs,e,5,2)]=P[:,14] #b41
    C[:,idx(gs,e,6,2)]=P[:,11] #b51

    C[:,idx(gs,e,1,3)]=P[:,7] #b02
    C[:,idx(gs,e,2,3)]=P[:,8] #b12
    C[:,idx(gs,e,3,3)]=P[:,9] #b22
    C[:,idx(gs,e,4,3)]=P[:,18] #b32
    C[:,idx(gs,e,5,3)]=P[:,15] #b42
    C[:,idx(gs,e,6,3)]=P[:,12] #b52

    C[:,idx(gs,e,1,4)]=P[:,30] #b03
    C[:,idx(gs,e,2,4)]=P[:,33] #b13
    C[:,idx(gs,e,3,4)]=P[:,36] #b23
    C[:,idx(gs,e,4,4)]=P[:,27] #b33
    C[:,idx(gs,e,5,4)]=P[:,26] #b43
    C[:,idx(gs,e,6,4)]=P[:,25] #b53

    C[:,idx(gs,e,1,5)]=P[:,29] #b04
    C[:,idx(gs,e,2,5)]=P[:,32] #b14
    C[:,idx(gs,e,3,5)]=P[:,35] #b24
    C[:,idx(gs,e,4,5)]=P[:,24] #b34
    C[:,idx(gs,e,5,5)]=P[:,23] #b44
    C[:,idx(gs,e,6,5)]=P[:,22] #b54

    C[:,idx(gs,e,1,6)]=P[:,28] #b05
    C[:,idx(gs,e,2,6)]=P[:,31] #b15
    C[:,idx(gs,e,3,6)]=P[:,34] #b25
    C[:,idx(gs,e,4,6)]=P[:,21] #b35
    C[:,idx(gs,e,5,6)]=P[:,20] #b45
    C[:,idx(gs,e,6,6)]=P[:,19] #b55

else

    D=fill(Vector{Float64}(),4,4);
    D[1,1]=P[:,1];
    D[1,2]=P[:,2];
    D[1,3]=P[:,7];
    D[1,4]=P[:,5];
    D[2,1]=P[:,3];
    D[2,2]=P[:,4];
    D[2,3]=P[:,8];
    D[2,4]=P[:,6];
    D[3,1]=P[:,14];
    D[3,2]=P[:,16];
    D[3,3]=P[:,12];
    D[3,4]=P[:,11];
    D[4,1]=P[:,13];
    D[4,2]=P[:,15];
    D[4,3]=P[:,10];
    D[4,4]=P[:,9];

    Q1=degree_elevate(D);
    Q2=degree_elevate(Q1);


    C[:,idx(gs,e,1,1)]=Q2[1] #b00
    C[:,idx(gs,e,1,2)]=Q2[2] #b10
    C[:,idx(gs,e,1,3)]=Q2[3] #b20
    C[:,idx(gs,e,1,4)]=Q2[4] #b30
    C[:,idx(gs,e,1,5)]=Q2[5] #b40
    C[:,idx(gs,e,1,6)]=Q2[6] #b50

    C[:,idx(gs,e,2,1)]=Q2[7] #b01
    C[:,idx(gs,e,2,2)]=Q2[8] #b11
    C[:,idx(gs,e,2,3)]=Q2[9] #b21
    C[:,idx(gs,e,2,4)]=Q2[10] #b31
    C[:,idx(gs,e,2,5)]=Q2[11] #b41
    C[:,idx(gs,e,2,6)]=Q2[12] #b51

    C[:,idx(gs,e,3,1)]=Q2[13] #b02
    C[:,idx(gs,e,3,2)]=Q2[14] #b12
    C[:,idx(gs,e,3,3)]=Q2[15] #b22
    C[:,idx(gs,e,3,4)]=Q2[16] #b32
    C[:,idx(gs,e,3,5)]=Q2[17] #b42
    C[:,idx(gs,e,3,6)]=Q2[18] #b52

    C[:,idx(gs,e,4,1)]=Q2[19] #b03
    C[:,idx(gs,e,4,2)]=Q2[20] #b13
    C[:,idx(gs,e,4,3)]=Q2[21] #b23
    C[:,idx(gs,e,4,4)]=Q2[22] #b33
    C[:,idx(gs,e,4,5)]=Q2[23] #b43
    C[:,idx(gs,e,4,6)]=Q2[24] #b53

    C[:,idx(gs,e,5,1)]=Q2[25] #b04
    C[:,idx(gs,e,5,2)]=Q2[26] #b14
    C[:,idx(gs,e,5,3)]=Q2[27] #b24
    C[:,idx(gs,e,5,4)]=Q2[28] #b34
    C[:,idx(gs,e,5,5)]=Q2[29] #b44
    C[:,idx(gs,e,5,6)]=Q2[30] #b54

    C[:,idx(gs,e,6,1)]=Q2[31] #b05
    C[:,idx(gs,e,6,2)]=Q2[32] #b15
    C[:,idx(gs,e,6,3)]=Q2[33] #b25
    C[:,idx(gs,e,6,4)]=Q2[34] #b35
    C[:,idx(gs,e,6,5)]=Q2[35] #b45
    C[:,idx(gs,e,6,6)]=Q2[36] #b55

end
end

#=

    C[:,idx(gs,e,1,1)]=P[:,1] #b00
    C[:,idx(gs,e,2,1)]=P[:,2] #b10
    C[:,idx(gs,e,3,1)]=P[:,7] #b20
    C[:,idx(gs,e,4,1)]=P[:,5] #b30

    C[:,idx(gs,e,1,2)]=P[:,3] #b01
    C[:,idx(gs,e,2,2)]=P[:,4] #b11
    C[:,idx(gs,e,3,2)]=P[:,8] #b21
    C[:,idx(gs,e,4,2)]=P[:,6] #b31

    C[:,idx(gs,e,1,3)]=P[:,14] #b02
    C[:,idx(gs,e,2,3)]=P[:,16] #b12
    C[:,idx(gs,e,3,3)]=P[:,12] #b22
    C[:,idx(gs,e,4,3)]=P[:,11] #b32

    C[:,idx(gs,e,1,4)]=P[:,13] #b03
    C[:,idx(gs,e,2,4)]=P[:,15] #b13
    C[:,idx(gs,e,3,4)]=P[:,10] #b23
    C[:,idx(gs,e,4,4)]=P[:,9] #b33

end
end




    m=gs.mesh;

    #I start ordering from e

    C[:,idx(gs,e,1,1)]=P[:,1] #b00
    C[:,idx(gs,e,1,2)]=P[:,2] #b10
    C[:,idx(gs,e,1,3)]=P[:,3] #b20
    C[:,idx(gs,e,2,1)]=P[:,4] #b01
    C[:,idx(gs,e,2,2)]=P[:,5] #b11
    C[:,idx(gs,e,2,3)]=P[:,6] #b21
    C[:,idx(gs,e,3,1)]=P[:,7] #b02
    C[:,idx(gs,e,3,2)]=P[:,8] #b12
    C[:,idx(gs,e,3,3)]=P[:,9] #b22

    e=next(m,e);

    C[:,idx(gs,e,1,1)]=P[:,10] #b50
    C[:,idx(gs,e,1,2)]=P[:,11] #b51
    C[:,idx(gs,e,1,3)]=P[:,12] #b52
    C[:,idx(gs,e,2,1)]=P[:,13] #b40
    C[:,idx(gs,e,2,2)]=P[:,14] #b41
    C[:,idx(gs,e,2,3)]=P[:,15] #b42
    C[:,idx(gs,e,3,1)]=P[:,16] #b30
    C[:,idx(gs,e,3,2)]=P[:,17] #b31
    C[:,idx(gs,e,3,3)]=P[:,18] #b32

    e=next(m,e);

    C[:,idx(gs,e,1,1)]=P[:,19] #b55
    C[:,idx(gs,e,1,2)]=P[:,20] #b45
    C[:,idx(gs,e,1,3)]=P[:,21] #b35
    C[:,idx(gs,e,2,1)]=P[:,22] #b54
    C[:,idx(gs,e,2,2)]=P[:,23] #b44
    C[:,idx(gs,e,2,3)]=P[:,24] #b34
    C[:,idx(gs,e,3,1)]=P[:,25] #b53
    C[:,idx(gs,e,3,2)]=P[:,26] #b43
    C[:,idx(gs,e,3,3)]=P[:,27] #b33

    e=next(m,e);

    C[:,idx(gs,e,1,1)]=P[:,28] #b05
    C[:,idx(gs,e,1,2)]=P[:,29] #b04
    C[:,idx(gs,e,1,3)]=P[:,30] #b03
    C[:,idx(gs,e,2,1)]=P[:,31] #b15
    C[:,idx(gs,e,2,2)]=P[:,32] #b14
    C[:,idx(gs,e,2,3)]=P[:,33] #b13
    C[:,idx(gs,e,3,1)]=P[:,34] #b25
    C[:,idx(gs,e,3,2)]=P[:,35] #b24
    C[:,idx(gs,e,3,3)]=P[:,36] #b23



end
=#
