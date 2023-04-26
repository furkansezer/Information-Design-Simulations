n=5; %number of agents

beta_f=norminv(0.99); %confidence levels
beta_l=norminv(0.99);

% delta=0.1:0.2:0.9;

% h_variance=0.1:0.1:1;
h_variance=1;

 
% app_factor=0.1;

solution=zeros(length(h_variance),2*n,2*n);

opt_value=zeros(length(h_variance),1);
S=zeros(2*n,2*n);
R=zeros(2*n,2*n);

coefficient_variance=zeros(2*n,2*n);

cov_theta=0.5*ones(n,n);

for i=1:n
    cov_theta(i,i)=5;
end 

    
Q=-1*ones(n);
for i=1:n
    Q(i,i)=20;
end 

no_info_solution=[zeros(n,n),zeros(n,n);zeros(n,n),Q];
full_info_solution=[Q\cov_theta/inv(Q).',Q\cov_theta;cov_theta/inv(Q).',cov_theta];

% V_ss=zeros(2*n,2*n); %agreement
%     for i=1:n
%         for j=1:n
%             V_ss(i,j)= 1/n;
%         end 
%             V_ss(i,i)=-1 + 1/n;
%     end

V_sw=[-Q eye(n);eye(n) zeros(n)];    %social welfare
    
for g=1:length(h_variance)
    overall_coefficient_variance=[h_variance(g)*ones(n) zeros(n);zeros(n) zeros(n)];
    
%     for w=1:5

     

cvx_begin quiet
variable X(2*n,2*n) symmetric semidefinite
variable t

maximize (t)
subject to

% trace((V_ss*X.'))>=t; %aggreement 

trace((V_sw*X.'))-beta_f*norm(trace(overall_coefficient_variance*X'))>=t; %social welfare


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
        coefficient_variance(k,i)=h_variance(g)/2;
        coefficient_variance(i,k)=h_variance(g)/2;
                
    end
    R(k,k)=Q(k,k);
    coefficient_variance(k,k)=h_variance(g);
    R(k,n+k)=-0.5;
    R(n+k,k)=-0.5;
    
    
    beta_l*norm(trace(coefficient_variance*X')) +trace(R*X') <=0; % second constraint group   
    

    R=zeros(2*n,2*n);
    coefficient_variance=zeros(2*n,2*n);  
   
end


% X>=0; %third constraint


cvx_end 

opt_value(g)=cvx_optval;
solution(g,:,:)=X;

%     end
end

% contourf(0.95:0.01:0.99,h_variance,opt_value)
 
pd = makedist('Normal',0,h_variance);
monte_h =random(pd,n,n,100);
monte_h=monte_h+Q;
%     