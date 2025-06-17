clc

figure(10)
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_xx/1000000,sigma_xx/1000000)
title('sxx')
axis equal
xlim([-0.04 0.04])
ylim([-0.04 0.04])
view(2)
shading interp
c=colorbar
title(c,'MPa')
caxis([-3e1 +3e1])
for i=1:nh
   plot3([nodes(5+(i-1)*n_ell,1);nodes(4+i*n_ell/2,1)],[nodes(5+(i-1)*n_ell,2);nodes(4+i*n_ell/2,2)], ...
       0.01*ones(2,1),'k','LineWidth',1)
end

figure(11)
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_yy/1000000,sigma_yy/1000000)
title('syy')
axis equal
view(2)
xlim([-0.04 0.04])
ylim([-0.04 0.04])
shading interp
c=colorbar
title(c,'MPa')
caxis([-3e1 +3e1])
for i=1:nh
   plot3([nodes(5+(i-1)*n_ell,1);nodes(4+i*n_ell/2,1)],[nodes(5+(i-1)*n_ell,2);nodes(4+i*n_ell/2,2)], ...
       0.01*ones(2,1),'k','LineWidth',1)
end

figure(12)
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_xy/1000000,sigma_xy/1000000)
title('sxy')
axis equal
view(2)
xlim([-0.04 0.04])
ylim([-0.04 0.04])
shading interp
c=colorbar
title(c,'MPa')
caxis([1e1 +3.8e1])
for i=1:nh
   plot3([nodes(5+(i-1)*n_ell,1);nodes(4+i*n_ell/2,1)],[nodes(5+(i-1)*n_ell,2);nodes(4+i*n_ell/2,2)], ...
       0.01*ones(2,1),'k','LineWidth',1)
end

figure(13)
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_1/1000000,sigma_1/1000000)
title('s1')
axis equal
view(2)
xlim([-0.04 0.04])
ylim([-0.04 0.04])
shading interp
c=colorbar
title(c,'MPa')
caxis([-0e1 +4.8e1])
for i=1:nh
   plot3([nodes(5+(i-1)*n_ell,1);nodes(4+i*n_ell/2,1)],[nodes(5+(i-1)*n_ell,2);nodes(4+i*n_ell/2,2)], ...
       0.01*ones(2,1),'k','LineWidth',1)
end

figure(14)
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_2/1000000,sigma_2/1000000)
title('s2')
axis equal
view(2)
xlim([-0.04 0.04])
ylim([-0.04 0.04])
shading interp
c=colorbar
title(c,'MPa')
caxis([-4.8e1 0e1])
for i=1:nh
   plot3([nodes(5+(i-1)*n_ell,1);nodes(4+i*n_ell/2,1)],[nodes(5+(i-1)*n_ell,2);nodes(4+i*n_ell/2,2)], ...
       0.01*ones(2,1),'k','LineWidth',1)
end

figure(15)
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*(sigma_1+sigma_2)/2/1000000,(sigma_1+sigma_2)/2/1000000)
title('sm')
axis equal
view(2)
xlim([-0.04 0.04])
ylim([-0.04 0.04])
shading interp
c=colorbar
title(c,'MPa')
caxis([-2.4e1 +2.4e1])
colormap('redblue')
for i=1:nh
   plot3([nodes(5+(i-1)*n_ell,1);nodes(4+i*n_ell/2,1)],[nodes(5+(i-1)*n_ell,2);nodes(4+i*n_ell/2,2)], ...
       0.01*ones(2,1),'k','LineWidth',1)
end




figure(16)
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*(sigma_1-sigma_2)/1000000,(sigma_1-sigma_2)/1000000)
title('sd')
axis equal
view(2)
xlim([-0.04 0.04])
ylim([-0.04 0.04])
shading interp
c=colorbar
title(c,'MPa')
caxis([2.4e1 +7.2e1])
for i=1:nh
   plot3([nodes(5+(i-1)*n_ell,1);nodes(4+i*n_ell/2,1)],[nodes(5+(i-1)*n_ell,2);nodes(4+i*n_ell/2,2)], ...
       0.01*ones(2,1),'k','LineWidth',1)
end



