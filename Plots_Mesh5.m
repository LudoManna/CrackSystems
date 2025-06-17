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
trisurf(elems(:,1:3),nodes1(:,1),nodes1(:,2),sigma_xx,sigma_xx)
title('sxx')
axis equal
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
view(2)
shading interp
colorbar
caxis([-0.2e8 +0.2e8]) %notnormalized
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('x(m)')
ylabel('y(m)')




figure(n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),sigma_yy,sigma_yy)
title('syy')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
colorbar
caxis([-0.2e8 +0.2e8]) %notnormalized
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('x(m)')
ylabel('y(m)')




figure(2*n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),sigma_xy,sigma_xy)
title('sxy')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
colorbar
caxis([0.1e8 +0.3e8]) %notnormalized
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('x(m)')
ylabel('y(m)')






figure(3*n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_1,sigma_1)
title('s1')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
colorbar
caxis([0.1e8 +0.4e8]) %notnormalized
xlabel('x(m)')
ylabel('y(m)')
% set(gca,'xtick',[])
% set(gca,'ytick',[])




figure(4*n^2+i+n*(j-1))
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_xy,sigma_2)
title('s2')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
colorbar
caxis([-0.45e8 0e8]) %notnormalized
xlabel('x(m)')
ylabel('y(m)')

nret=100;
yzeros=linspace(lims(2),lims(1),nret)';
xzeros=lims(1)*ones(nret,1);
bs=yzeros+xzeros*tan(pi/12);
xones=ones(nret,1)*lims(2);
yones=-xones*tan(pi/12)+bs;
for iret=1:nret
    plot([xzeros(iret) xones(iret)],[yzeros(iret) yones(iret)],'r')
end
% quiver(nodes(:,1),nodes(:,2),cos(theta_sig1),sin(theta_sig1),1)
% set(gca,'xtick',[])
% set(gca,'ytick',[])

figure(10)
hold on
% plot([nodes(1:4,1);nodes(1,1)],[nodes(1:4,2);nodes(1,2)],'k')
% title('s2 max')
% axis equal
% for ic=1:n_crack
%     plot([nodes(5+n_ell*(ic-1):4+n_ell*ic,1);nodes(5+n_ell*(ic-1),1)], ...
%         [nodes(5+n_ell*(ic-1):4+n_ell*ic,2);nodes(5+n_ell*(ic-1),2)],'k')
% end
s2min=min(sigma_2);
% sxymax_i=find(sigma_xy==sxymax);
s2min_i=find(sigma_2<=s2min*2/3);
scatter(nodes(s2min_i,1),nodes(s2min_i,2),90,'k','filled')
scatter(nodes(s2min_i,1),nodes(s2min_i,2),50,'w','filled')
% xlabel('x(m)')
% ylabel('y(m)')

figure(10)
hold on
plot([nodes(1:4,1);nodes(1,1)],[nodes(1:4,2);nodes(1,2)],'k')
% title('sxy max')
axis equal
for ic=1:n_crack
    plot([nodes(5+n_ell*(ic-1):4+n_ell*ic,1);nodes(5+n_ell*(ic-1),1)], ...
        [nodes(5+n_ell*(ic-1):4+n_ell*ic,2);nodes(5+n_ell*(ic-1),2)],'k')
end
sxymax=max(sigma_xy);
% sxymax_i=find(sigma_xy==sxymax);
sxymax_i=find(sigma_xy>=sxymax*3/4);
scatter(nodes(sxymax_i,1),nodes(sxymax_i,2),25,'r','filled')
xlabel('x(m)')
ylabel('y(m)')
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])

% figure(11)
% hold on
% plot([nodes(1:4,1);nodes(1,1)],[nodes(1:4,2);nodes(1,2)],'k')
% % title('sxy max')
% axis equal
% for ic=1:n_crack
%     plot([nodes(5+n_ell*(ic-1):4+n_ell*ic,1);nodes(5+n_ell*(ic-1),1)], ...
%         [nodes(5+n_ell*(ic-1):4+n_ell*ic,2);nodes(5+n_ell*(ic-1),2)],'k')
% end
% scatter(nodes(5:n_ell/2:5+n_ell*(n_crack-1)+n_ell/2,1),nodes(5:n_ell/2:5+ ...
%     n_ell*(n_crack-1)+n_ell/2,2),25,'r','filled')
% xlabel('x(m)')
% ylabel('y(m)')
% xlim([lims(i) lims(i+1)])
% ylim([lims(j) lims(j+1)])
ALP=zeros(2*n_crack,1);
MS2=zeros(n_ell,n_crack);
MSxy=zeros(n_ell,n_crack);

for ialp=1:n_crack
    ALP(2*ialp-1,1)=alpha(ialp)*180/pi;
    ALP(2*ialp,1)=alpha(ialp)*180/pi;
    MS2(:,ialp)=-sigma_2(5+n_ell*(ialp-1):4+n_ell*ialp);
        MSxy(:,ialp)=sigma_xy(5+n_ell*(ialp-1):4+n_ell*ialp);

end
MS2=max(MS2);
MS2=MS2';
MSxy=max(MSxy);
MSxy=MSxy';

figure(12)
hold on

scatter(ALP,-sigma_2(5:n_ell/2:5+ ...
    n_ell*(n_crack-1)+n_ell/2),25,'r','filled')

scatter(alpha*180/pi,sigma_xy(5:n_ell:5+ ...
    n_ell*(n_crack-1)),25,'b','filled')
plot([min(ALP);max(ALP)],[-s2min;-s2min],'r')
plot([min(ALP);max(ALP)],[sxymax;sxymax],'b')

figure(13)
hold on

scatter(alpha*180/pi,MS2,25,'r','filled')

scatter(alpha*180/pi,MSxy,25,'b','filled')
plot([min(ALP);max(ALP)],[-17/24*s2min;-17/24*s2min],'r')
plot([min(ALP);max(ALP)],[17/24*sxymax;17/24*sxymax],'b')
% xlim([min(ALP) max(ALP)])
% ylim([lims(j) lims(j+1)])
pbaspect([1 1 1])


% 
figure(5*n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_1,(sigma_1+sigma_2)/2)
title('sm')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
colorbar
caxis([-0.2e8 +0.2e8]) %notnormalized
colormap('redblue')
xlabel('x(m)')
ylabel('y(m)')
% set(gca,'xtick',[])
% set(gca,'ytick',[])


% 
% 
% 
figure(6*n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*sigma_xy,sigma_1-sigma_2)
title('sd')
axis equal
view(2)
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
shading interp
colorbar
caxis([0.4e8 +0.7e8]) %notnormalized
xlabel('x(m)')
ylabel('y(m)')
% set(gca,'xtick',[])
% set(gca,'ytick',[])


figure(7*n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*uhx,uhx)
title('ux')
axis equal
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
view(2)
shading interp
colorbar
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('x(m)')
ylabel('y(m)')



figure(8*n^2+i+n*(j-1))
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*uhy,uhy)
title('uy')
axis equal
xlim([lims(i) lims(i+1)])
ylim([lims(j) lims(j+1)])
view(2)
shading interp
colorbar
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('x(m)')
ylabel('y(m)')



end
end





