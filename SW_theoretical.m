n=4;
% opt_value=zeros(1,496);
% var_a1=zeros(496,1);
% cov_a1_a2=zeros(496,1);
% f=zeros(1,496);
% var_t1=zeros(496,1);
% cov_t1_t2=zeros(496,1);
% t=0:0.002:0.99;
t=0.2;
z=1;

cov_theta=ones(n,n);
    for i=1:n
        cov_theta(i,i)=2;
    end 
    
% for z=1:496
Q=ones(n,n);
Q=t(z)*Q;   


    for i=1:n
        Q(i,i)=1; 
    end
f(z)=trace(Q\cov_theta');
S=zeros(2*n,2*n);
R=zeros(2*n,2*n);

V_sw=[-Q,eye(4);eye(4),zeros(4)];

% Sum of Squares objective


cvx_begin sdp quiet
variable X(2*n,2*n) symmetric
dual variable mu{n,n}

maximize (trace(V_sw*X'))

for i=1:n
    for j=1:i
        if(i==j)
            S(n+i,n+j)=1;
        mu{i,j}: S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %first constraint
            S(n+i,n+j)=0;
            
        end
        if (j<i)
            S(n+i,n+j)=1;
        mu{i,j}: S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %first constraint
       
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

var_a1(z)=X(1,1);
cov_a1_a2(z)=X(1,2);


%Calculation of cov(t, theta)
M=zeros(n);
for i=1:n
    for j=1:n
        M(i,j)=Q(i,j)/Q(i,i);
    end
end
for i=1:n
   M(i,i)=-1; 
end
N=diag(Q);
b=linsolve(M,N);


cov_t_theta=X;
for i=1:n %calculation of cov(t_i, theta_i)
    for j=n+1:2*n
     cov_t_theta(i,j)= X(i,j)/b(i);  
     cov_t_theta(j,i)= X(j,i)/b(j-n);
    end
end
for i=1:n %calculation of cov(t_i, t_j)
    for j=1:n
        cov_t_theta(i,j)=X(i,j)/(b(i)*b(j));
    end
end
var_t1(z)=cov_t_theta(1,1);
cov_t1_t2(z)=cov_t_theta(1,2);



opt_value(z)=cvx_optval;
%  end
% plot(t,var_a1);
% plot(t,var_t1);
% plot(t,cov_a1_a2);
% plot(t,cov_t1_t2);
% corr_a=cov_a1_a2./var_a1;
% corr_t=cov_t1_t2./var_t1;
% plot(t,corr_a)
% ylim([0 2])
% plot(t,corr_t) 
% ylim([0 2])
