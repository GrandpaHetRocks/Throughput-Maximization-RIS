function[hij_array,theta,Hri_array,G_array,vij_array,Vri_array]=channel_generation2(k,n)

    Kg=3;
    Kh=3;
    Kv=3;
    dh=5;
    dv=5;
    dg=5;
    alpha=2.2;
    alphad=3;
    dd=10;
    beta_g=10^-3*dg^(-alpha);
    beta_h=10^-3*dh^(-alpha);
    beta_v=10^-3*dv^(-alpha);
    beta_direct=10^-3*dd^(-alphad);

    
    a=randn(k,k);
    b=randn(k,k);
    hij_array=(sqrt(0*beta_direct/(1+0)))+(sqrt(beta_direct/((0+1)*2)))*complex(a,b);

    a=randn(k,k);
    b=randn(k,k);
    vij_array=(sqrt(0*beta_direct/(1+0)))+(sqrt(beta_direct/((0+1)*2)))*complex(a,b);

    a=randn(1,n);
    b=randn(1,n);
    theta=complex(a,b)';
    theta = bsxfun(@rdivide, theta, sqrt(sum(theta.*conj(theta), 2)));

    Hri_array=[];
    for i=1:k
        a=randn(n,1);
        b=randn(n,1);
        Hri=(sqrt(Kh*beta_h/(1+Kh)))+(sqrt(beta_h/((Kh+1)*2)))*complex(a,b);
        Hri_array=[Hri_array Hri];
    end
    
    Vri_array=[];
    for i=1:k
        a=randn(n,1);
        b=randn(n,1);
        Vri=(sqrt(Kv*beta_v/(1+Kv)))+(sqrt(beta_v/((Kv+1)*2)))*complex(a,b);
        Vri_array=[Vri_array Vri];
    end

    G_array=[];
    for i=1:k
        a=randn(n,1);
        b=randn(n,1);
        G=(sqrt(Kg*beta_g/(1+Kg)))+(sqrt(beta_g/((Kg+1)*2)))*complex(a,b);
        G_array=[G_array G];
    end


%     v_array=[];
%     for i=1:k
%         a=randn(n,1);
%         b=randn(n,1);
%         V=complex(a,b);
%         v_array=[v_array V];
%     end
end