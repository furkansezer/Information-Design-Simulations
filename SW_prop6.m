n=4;

a=0.5:0.01:0.99;

cov_theta=ones(n,n);
Q=zeros(30,n,n);
f=zeros(30,50);

opt_value=zeros(30,50);
var_a1=zeros(30,50);
cov_a1_a2=zeros(30,50);

% opt_value=zeros(50,1);
% var_a1=zeros(50,1);
% cov_a1_a2=zeros(50,1);
% f=zeros(50,1);
 
for m=1:30    
    
    r=rand(n);
%     Q(m,:,:)=r+r';
    Q(m,:,:)=r;
    for i=1:n
        Q(m,i,i)=4;
    end    

     for z=1:50
    
        for i=1:n
        cov_theta(i,i)=1/a(z);
        end 
         Q_temp=reshape(Q(m,:,:),[4,4]);
        f(m,z)=trace(Q_temp\cov_theta'); %full info objective value
    
S=zeros(2*n,2*n);
R=zeros(2*n,2*n);

V_sw=[-Q_temp,eye(4);eye(4),zeros(4)];

% Sum of Squares objective


cvx_begin sdp quiet
variable X(2*n,2*n) symmetric
% dual variable mu{n,n}

maximize (trace(V_sw*X'))

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
        R(k,i)=Q_temp(k,i)/2;
        R(i,k)=Q_temp(k,i)/2;
    end
    R(k,k)=Q_temp(k,k);
    R(k,n+k)=-0.5;
    R(n+k,k)=-0.5;
    trace(R*X')==0; % second constraint
    R=zeros(2*n,2*n);    
    
end

X>=0; %third constraint

cvx_end 

var_a1(m,z)=X(1,1);
cov_a1_a2(m,z)=X(1,2);

opt_value(m,z)=cvx_optval;

    end
end

avg_opt_value=mean(opt_value,1);
avg_f=mean(f,1);
% plot(a,avg_opt_value,a,avg_f)
% plot(a,opt_value,a,f)