n=5;
eps2=0.0001;
eps1=0.03:0.01:0.12;
rho=0.7:0.3:3.4;

% eps2=0.01;
% eps1=0.01:0.10:0.91;
% rho=0.3:0.15:1.65;

app_factor=0.1;

solution=zeros(10,10,2*n,2*n);

opt_value=zeros(10,10);
S=zeros(2*n,2*n);
R=zeros(2*n,2*n);

shift=zeros(2*n,2*n);

cov_theta=0.5*ones(n,n);

for i=1:n
    cov_theta(i,i)=5;
end 

    
Q=-1*ones(n);
for i=1:n
    Q(i,i)=5;
end 

no_info_solution=[zeros(n,n),zeros(n,n);zeros(n,n),cov_theta];
full_info_solution=[Q\cov_theta/inv(Q).',Q\cov_theta;cov_theta/inv(Q).',cov_theta];

V_ss=zeros(2*n,2*n); %agreement
    for i=1:n
        for j=1:n
            V_ss(i,j)= 1/n;
        end 
            V_ss(i,i)=-1 + 1/n;
    end

V_sw=[-Q eye(n);eye(n) zeros(n)];    %social welfare
    
for g=1:10
    overall_shift=[eps2*ones(n) zeros(n);zeros(n) zeros(n)];
    for i=1:n
        overall_shift(i,i)=eps1(g);
    end
    
    for w=1:10

     

cvx_begin quiet
variable X(2*n,2*n) symmetric semidefinite
variable t

maximize (t)
subject to

% trace((V_ss*X.'))>=t; %aggreement 

trace((V_sw*X.'))-((n^2)*rho(w)/(2*n-1))*norm(trace(overall_shift*X'))>=t; %social welfare


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
        shift(k,i)=eps2/2;
        shift(i,k)=eps2/2;
                
    end
    R(k,k)=Q(k,k);
    shift(k,k)=eps1(g);
    R(k,n+k)=-0.5;
    R(n+k,k)=-0.5;
    
    
    rho(w)*norm(trace(shift*X')) +trace(R*X') <=app_factor; % second constraint group   
    rho(w)*norm(trace(shift*X')) -trace(R*X') <=app_factor; % second constraint group
    

    R=zeros(2*n,2*n);
    shift=zeros(2*n,2*n);  
    
end


% X>=0; %third constraint


cvx_end 

opt_value(g,w)=cvx_optval;
solution(g,w,:,:)=X;

    end
 end

dist_no=zeros(10,10);
dist_full=zeros(10,10);
dist_full_no=sqrt(sum((full_info_solution-no_info_solution).^2,'all'));

for i=1:10
    for j=1:10
        dist_no(i,j)=sqrt(sum((squeeze(solution(i,j,:,:))-no_info_solution).^2,'all'));
        dist_full(i,j)=sqrt(sum((squeeze(solution(i,j,:,:))-full_info_solution).^2,'all'));
    end    
end    
 
% normalized_dist_no=dist_no./max(dist_no, [], 'all'); %divide by max element in dist no
% contourf(eps1,rho,normalized_dist_no)


total_normalized_dist_no=dist_no./dist_full_no;
% total_normalized_dist_full=dist_full./dist_full_no;
contourf(eps1,rho,total_normalized_dist_no) 
% contourf(eps1,rho,total_normalized_dist_full)
% contourf(eps1,rho,opt_value) 

    