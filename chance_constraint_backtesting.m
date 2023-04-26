
coefficient_variance=zeros(2*n,2*n);
R=zeros(2*n,2*n);
count=0;
g=1;
    k=1;
for j=1:100
    for i=1:n
        R(k,i)=squeeze(monte_h(k,i,j))/2;
        R(i,k)=squeeze(monte_h(k,i,j))/2;
        coefficient_variance(k,i)=h_variance(g)/2;
        coefficient_variance(i,k)=h_variance(g)/2;
                
    end
    R(k,k)=squeeze(monte_h(k,k,j));
    coefficient_variance(k,k)=h_variance(g);
    R(k,n+k)=-0.5;
    R(n+k,k)=-0.5;
    
    
    if(beta_l*norm(trace(coefficient_variance*X')) +trace(R*X') <=0.1) 
        count=count+1;
    end    
    
    R=zeros(2*n,2*n);
    coefficient_variance=zeros(2*n,2*n); 
end
