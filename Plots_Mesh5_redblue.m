clc
theta_sig1=atan(sigmavec_1(:,2)./sigmavec_1(:,1)); %radians
n=1;
lims=zeros(n+1,1);
lims(1)=-box_a;
for i=2:n+1
    lims(i)=lims(i-1)+2*box_a/n;
end

for j=1:n
for i=1:n
figure(i+n*(j-1))
trisurf(elems(:,1:3),nodes1(:,1),nodes1(:,2),sigma_xx,-sigma_xx/1000000)
% title('sxx')
axis equal
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
view(2)
shading interp
c=colorbar
title(c,'MPa')
colormap('redblue')
caxis([-0.2e2 +0.2e2]) %notnormalized
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('x(m)')
ylabel('y(m)')




figure(n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),sigma_yy,-sigma_yy/1000000)
% title('syy')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
c=colorbar
title(c,'MPa')
colormap('redblue')
caxis([-0.1e2 +0.1e2]) %notnormalized
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('x(m)')
ylabel('y(m)')




figure(2*n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),sigma_xy,sigma_xy/1000000)
% title('sxy')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
c=colorbar
title(c,'MPa')
colormap('redblue')
caxis([0.15e2 +0.3e2]) %notnormalized
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('x(m)')
ylabel('y(m)')






figure(3*n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_1,sigma_1/1000000)
% title('s1')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
c=colorbar
title(c,'MPa')
colormap('redblue')
caxis([-0.55e2 +0.55e2]) %notnormalized
xlabel('x(m)')
ylabel('y(m)')
% set(gca,'xtick',[])
% set(gca,'ytick',[])




figure(4*n^2+i+n*(j-1))
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_xy,sigma_2/1000000)
% title('s2')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
% colorbar
c=colorbar
title(c,'MPa')
colormap('redblue')
caxis([-0.55e2 0.55e2]) %notnormalized
xlabel('x(m)')
ylabel('y(m)')


figure(5*n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_1,(sigma_1+sigma_2)/2/1000000)
% title('sm')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
c=colorbar
title(c,'MPa')
colormap('redblue')
caxis([-0.2e2 +0.2e2]) %notnormalized
xlabel('x(m)')
ylabel('y(m)')
% set(gca,'xtick',[])
% set(gca,'ytick',[])


% 
% 
% 
figure(6*n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_xy,(sigma_1-sigma_2)/1000000)
% title('sd')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
c=colorbar
title(c,'MPa')
colormap('redblue')
caxis([0.25e2 +0.65e2]) %notnormalized
xlabel('x(m)')
ylabel('y(m)')
% set(gca,'xtick',[])
% set(gca,'ytick',[])




end
end





