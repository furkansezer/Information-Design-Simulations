%Lambda \in (0, 1) is the weighting parameter of SW with respect to SSD.
n=4;
lambda=0:0.01:1;

opt_value=zeros(1,101);
var_a1=zeros(1,101);
cov_a1_a2=zeros(1,101);
full_info_obj=zeros(1,101);

y=ones(4);    
S=zeros(2*n,2*n);
R=zeros(2*n,2*n);

cov_theta=ones(n,n);

    for i=1:n
        cov_theta(i,i)=4;
    end 

    
Q=-0.3*ones(4);

for i=1:n
    Q(i,i)=1;
end    

% Sum of Squares + \lambda*Social welfare objective
V_ss=zeros(2*n,2*n);
    for i=1:n
        for j=1:n
            V_ss(i,j)= 1/n;
        end 
            V_ss(i,i)=-1 + 1/n;
    end



V_sw=[-Q,eye(4);eye(4),zeros(4)];
    
for t=1:101
cvx_begin sdp quiet
variable X(2*n,2*n) symmetric

maximize (trace(((1-lambda(t))*V_sw +lambda(t)*V_ss)*X'))

for i=1:n
    for j=1:i
        if(i==j)
            S(n+i,n+j)=1;
            S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %first constraint
            S(n+i,n+j)=0;
            
        end
        if (j<i)
            S(n+i,n+j)=1;
           S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %first constraint
       
             S(n+i,n+j)=0;
        end
    end
end

for k=1:n
    for i=1:n
        R(k,i)=Q(k,i)/2;
        R(i,k)=Q(k,i)/2;
    end
    R(k,k)=Q(k,k);
    R(k,n+k)=-0.5;
    R(n+k,k)=-0.5;
    trace(R*X')==0; % second constraint
        
    R=zeros(2*n,2*n);
end

X>=0; %third constraint

cvx_end 
var_a1(t)=X(1,1);
cov_a1_a2(t)=X(1,2);

opt_value(t)=cvx_optval;

for i=1:4
    for j=1:4
        if (i==j)
            y(i,j)=lambda(t)/4-1;
        else
            y(i,j)=-0.05*lambda(t)+0.3;  
        end
    end
end    
v_q=inv(Q)'*(y+2*(1-lambda(t))*Q)*inv(Q); 

full_info_obj(t)=trace(v_q*cov_theta);
end

plot(lambda,opt_value,lambda,full_info_obj,lambda,zeros(1,101),'LineWidth',3)

