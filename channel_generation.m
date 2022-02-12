function[hij_array,theta,Hri_array,G_array]=channel_generation(k,n)
    a=randn(k,k);
    b=randn(k,k);
    hij_array=complex(a,b);
    a=randn(1,n);
    b=randn(1,n);
    theta=complex(a,b)';
    theta = bsxfun(@rdivide, theta, sqrt(sum(theta.*conj(theta), 2)));
    Hri_array=[];
    for i=1:k
        a=randn(n,1);
        b=randn(n,1);
        Hri=complex(a,b);
        Hri_array=[Hri_array Hri];
    end

    a=randn(1,k);
    b=randn(1,k);
    G_array=complex(a,b);
end