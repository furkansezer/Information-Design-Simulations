n=5; %number of agents
rng(0,'twister');
J=30; %number of samples

beta_f=0.95; %confidence levels
beta_l=0.95;

nu_f=((1-beta_f)*J)^(-1);
nu_l=((1-beta_l)*J)^(-1);

h_variance=0.1:0.1:1;

% app_factor=0.1;

solution=zeros(10,2*n,2*n);

opt_value=zeros(10,1);
S=zeros(2*n,2*n);
R=zeros(2*n,2*n);

coefficient_variance=zeros(2*n,2*n);

cov_theta=0.5*ones(n,n);

for i=1:n
    cov_theta(i,i)=5;
end 

    
% no_info_solution=[zeros(n,n),zeros(n,n);zeros(n,n),Q];
% full_info_solution=[Q\cov_theta/inv(Q).',Q\cov_theta;cov_theta/inv(Q).',cov_theta];

% V_ss=zeros(2*n,2*n); %agreement
%     for i=1:n
%         for j=1:n
%             V_ss(i,j)= 1/n;
%         end 
%             V_ss(i,i)=-1 + 1/n;
%     end


    
for g=1:10
    
    Q=-1*ones(n);  %Q mean
        for i=1:n
            Q(i,i)=20;
        end
         
    Q=h_variance(g).*randn(n,n,J)+Q;
    

V_sw=zeros(2*n,2*n,J);
for b=1:J
    
V_sw(:,:,b)=[squeeze(-Q(:,:,b)) eye(n);eye(n) zeros(n)];    %social welfare    %need to solve this

end



     

cvx_begin quiet
variable X(2*n,2*n) symmetric semidefinite
variable t
variable rho_f
variable rho_l(n)
variable z_f(J) nonnegative
variable z_l(n,J) nonnegative


maximize (t)
subject to

for b=1:J
    
rho_f+ nu_f*sum(z_f)<=0; %social welfare obj

z_f(b)>=t-trace(squeeze(V_sw(:,:,b))*X')-rho_f; %social welfare obj


for k=1:n
    for i=1:n
        R(k,i)=Q(k,i,b)/2;
        R(i,k)=Q(i,k,b)/2;

                
    end
    R(k,k)=Q(k,k,b);
    R(k,n+k)=-0.5;
    R(n+k,k)=-0.5;
    
    
    rho_l(k)+nu_l*sum(z_l(k,b))<=0;            % equilibrium constraint  
    z_l(k,b)>=trace(R*X')-rho_l(k);            % equilibrium constraint   
    

    R=zeros(2*n,2*n);
   
end

end


for i=1:n
    for j=1:i
        if(i==j)
            S(n+i,n+j)=1;
            S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %assignment constraint
            S(n+i,n+j)=0;
            
        end
        if (j<i)
            S(n+i,n+j)=1;
           S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %assignment constraint
       
             S(n+i,n+j)=0;
        end
    end
end




cvx_end 

opt_value(g)=cvx_optval;
solution(g,:,:)=X;

    
end

% contourf(0.95:0.01:0.99,h_variance,opt_value)
 
% plot(h_variance,opt_value,'LineWidth',3)
    