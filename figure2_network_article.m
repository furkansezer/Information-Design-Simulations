n=3:20;
z=1./(n-1);

u=(2*(n-1)+sqrt(n.^2-2*n+4))./(n.*(n-2));
l=(2*(n-1)-sqrt(n.^2-2*n+4))./(n.*(n-2));

plot(n,z,n,u,n,l)