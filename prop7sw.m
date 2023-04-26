n=4;

opt_value=zeros(4,17);
% var_a1=zeros(17,1);
% cov_a1_a2=zeros(17,1);
f=zeros(4,17);

h=0.2:0.2:0.8;
% h=-0.3:0.1:-0.1;

v=0.4:0.005:0.48;

for p=1:4
    
Q=h(p)*ones(n,n);
 
for i=1:n
        Q(i,i)=1; 
end    
for z=1:17
    
    cov_theta=0.2*ones(n,n);
    for i=1:n
        cov_theta(i,i)=v(z);
    end 

f(p,z)=trace(Q\cov_theta');
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

% var_a1(z)=X(1,1);
% cov_a1_a2(z)=X(1,2);

opt_value(p,z)=cvx_optval;
end
end

% plot(v,opt_value,v,f,v,zeros(1,17),'LineWidth',3)

hold on
plot(v,opt_value(1,:),'lineWidth',3)
hold on
plot(v,opt_value(2,:),'lineWidth',3)
hold on
plot(v,opt_value(3,:),'lineWidth',3)
hold on
plot(v,opt_value(4,:),'lineWidth',3)
% hold on
% plot(v,opt_value(5,:),'lineWidth',3)
hold on
plot(v,zeros(1,17),'LineWidth',3)
