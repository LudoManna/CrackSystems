% clc
% 
% figure(10)
% hold on
% trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_xx,sigma_xx/1000000)
% % title('\sigma_{xx}')
% axis equal
% xlim([-0.07 0.07])
% ylim([-0.01 0.01])
% view(2)
% shading interp
% colorbar
% % c=colorbar;
% % title(c,'MPa')
% caxis([-4e1 +4e1])
% for i=1:nh
%    plot(nodes(5+(i-1)*n_ell:4+i*n_ell,1),nodes(5+(i-1)*n_ell:4+i*n_ell,2),'k')
% end
% % set(gca,'xtick',[])
% 
% 
% 
% figure(11)
% hold on
% trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_yy,sigma_yy/1000000)
% % title('\sigma_{yy}')
% axis equal
% view(2)
% xlim([-0.07 0.07])
% ylim([-0.01 0.01])
% shading interp
% colorbar
% % c=colorbar;
% % title(c,'MPa')
% caxis([-40 +40])
% for i=1:nh
%    plot(nodes(5+(i-1)*n_ell:4+i*n_ell,1),nodes(5+(i-1)*n_ell:4+i*n_ell,2),'k')
% end
% % set(gca,'xtick',[])
% 
% 

% figure(12)
% hold on
% trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_xy,sigma_xy/1000000)
% % title('\sigma_{xy}')
% axis equal
% view(2)
% xlim([-0.07 0.07])
% ylim([-0.01 0.01])
% shading interp
% colorbar
% % c=colorbar;
% % title(c,'MPa')
% caxis([12 36])
% for i=1:nh
%    plot(nodes(5+(i-1)*n_ell:4+i*n_ell,1),nodes(5+(i-1)*n_ell:4+i*n_ell,2),'k')
% end
% % set(gca,'xtick',[])

% 
% figure(13)
% hold on
% trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_1,sigma_1/1000000)
% % title('\sigma_{1}')
% axis equal
% view(2)
% xlim([-0.07 0.07])
% ylim([-0.01 0.01])
% shading interp
% colorbar
% % c=colorbar;
% % title(c,'MPa')
% caxis([0 48])
% for i=1:nh
%    plot(nodes(5+(i-1)*n_ell:4+i*n_ell,1),nodes(5+(i-1)*n_ell:4+i*n_ell,2),'k')
% end
% % set(gca,'xtick',[])
Fi2=33*pi/180;
nodes2=nodes*[cos(Fi2) -sin(Fi2); sin(Fi2) cos(Fi2)];
sigmavec_1_2=sigmavec_1*[cos(Fi2) -sin(Fi2); sin(Fi2) cos(Fi2)];
figure(14)
hold on
% trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_2,sigma_2/1000000)
trisurf(elems(:,1:3),nodes2(:,1),nodes2(:,2),0*sigma_2,sigma_2/1000000)
% scatter(nodes2(:,1),nodes2(:,2),0.1,'k','filled')
% title('\sigma_{2}')
axis equal
view(2)
xlim([-0.12 0.12])
ylim([-0.06 0.06])
% xlim([-2 0])
% ylim([-0.5 0.5])
shading interp
colorbar
% c=colorbar;
% title(c,'MPa')
caxis([-48 0])
% quiver(nodes(5+n_ell*nh:end,1),nodes(5+n_ell*nh:end,2),sigmavec_1(5+n_ell*nh:end,1),sigmavec_1(5+n_ell*nh:end,2),0.01,'w','ShowArrowHead','off','linewidth',0.01)
for i=1:nh
   plot(nodes2(5+(i-1)*n_ell:4+i*n_ell,1),nodes2(5+(i-1)*n_ell:4+i*n_ell,2),'k')
end
set(gca,'xtick',[])
set(gca,'ytick',[])


% quiver(nodes2(:,1),nodes2(:,2),sigmavec_1_2(:,1),sigmavec_1_2(:,2),0.1,'w','ShowArrowHead','off')
% 
% % figure(15)
% % hold on
% % trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*(sigma_1+sigma_2)/2,(sigma_1+sigma_2)/2/1000000)
% % % title('sm')
% % axis equal
% % view(2)
% % xlim([-0.07 0.07])
% % ylim([-0.01 0.01])
% % shading interp
% % colorbar
% % caxis([-24 24])
% % colormap('redblue')
% % for i=1:nh
% %    plot(nodes(5+(i-1)*n_ell:4+i*n_ell,1),nodes(5+(i-1)*n_ell:4+i*n_ell,2),'k')
% % end
% 
% 
% 
% 
figure(16)
hold on
trisurf(elems(:,1:3),nodes2(:,1),nodes2(:,2),0*(sigma_1-sigma_2),(sigma_1-sigma_2)/1000000)
% title('\sigma_{d}')
axis equal
view(2)
xlim([-0.12 0.12])
ylim([-0.06 0.06])
shading interp
colorbar
% c=colorbar;
% title(c,'MPa')
caxis([-24 24])
caxis([20 76])
for i=1:nh
   plot(nodes2(5+(i-1)*n_ell:4+i*n_ell,1),nodes2(5+(i-1)*n_ell:4+i*n_ell,2),'k')
end
set(gca,'xtick',[])
set(gca,'ytick',[])

% 
% % quiver(nodes(5+n_ell*nh:end,1),nodes(5+n_ell*nh:end,2),sigmavec_1(5+n_ell*nh:end,1),sigmavec_1(5+n_ell*nh:end,2),0.05,'r','ShowArrowHead','off','linewidth',0.01)
% 
% 
% 
% 
% 
