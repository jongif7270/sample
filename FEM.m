function [A] = FEM(xl,xr,yl,yr,Mx,My,N)
[c4n,n4e,ind4e,inddb] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,N);
[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V);
b = zeros(size(c4n,1),1);
u = b;
f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
u_D=@(x) x(:,1)*0;

I = eye((N+1)*(N+2)/2);
xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

Kr=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2,size(n4e,1));

for j=1:size(n4e,1)
M=J(j)*I/(V*V');
Srr=J(j)*(V\Dr)'*(V\Dr);
Srs=J(j)*(V\Dr)'*(V\Ds);
Ssr=J(j)*(V\Ds)'*(V\Dr);
Sss=J(j)*(V\Ds)'*(V\Ds);

K=J(j)*((rx(j)^2+ry(j)^2)*Srr+(rx(j)*sx(j)+ry(j)*sy(j))*(Srs+Ssr)+(sx(j)^2+sy(j)^2)*Sss);
Kr(:,:,j)=K;

b(ind4e(j,:)) = b(ind4e(j,:)) + J(j)*M*f(c4n(ind4e(j,:),:));
end

ind=ind4e';
Ir=repmat(ind4e,1,(N+1)*(N+2)/2)';
Jr=(repmat(ind(:),1,(N+1)*(N+2)/2))';
A=sparse(Ir(:),Jr(:),Kr(:));

fns = setdiff(1:size(c4n,1), inddb);
u(inddb) = u_D(c4n(inddb,:));
u(fns) = A(fns,fns)\b(fns);
plot3(c4n(:,1),c4n(:,2),u)