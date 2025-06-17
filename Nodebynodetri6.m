function [sigma_xx,sigma_yy,sigma_xy,sigma,sigma_1,sigma_2,sigmavec_1,sigmavec_2]=Nodebynodetri6(nnod,elems,Sigma_xx,Sigma_yy,Sigma_xy)
%occhio a sigma, sigma_xx ,sigma_yy e sigma_xy
cornernodes=zeros(nnod,1);
sigma_xx=cornernodes;
sigma_yy=sigma_xx;
sigma_xy=sigma_xx;
sigma=zeros(2,2,nnod);
sigmavec_1=zeros(nnod,2);
sigmavec_2=zeros(nnod,2);
sigma_1=zeros(nnod,1);
sigma_2=zeros(nnod,1);
tic
for i=1:nnod
    element=find(elems(:,1)==i);
    element=[element,ones(size(element,1),1)];
    element=[element;[find(elems(:,2)==i),2*ones(size(find(elems(:,2)==i)))]];
    element=[element;[find(elems(:,3)==i),3*ones(size(find(elems(:,3)==i)))]];
    
    element1=[find(elems(:,4)==i),4*ones(size(find(elems(:,4)==i)))];
    element1=[element1;[find(elems(:,5)==i),5*ones(size(find(elems(:,5)==i)))]];
    element1=[element1;[find(elems(:,6)==i),6*ones(size(find(elems(:,6)==i)))]];
    
    neli=size(element,1);
    neli1=size(element1,1);
    if neli==0
        for ieli1=1:neli1
            if element1(ieli1,2)==4
            sigma_xx(i)=sigma_xx(i)+0.5*(sigma_xx(elems(element1(ieli1,1),2))+sigma_xx(elems(element1(ieli1,1),3)));
            sigma_yy(i)=sigma_yy(i)+0.5*(sigma_yy(elems(element1(ieli1,1),2))+sigma_yy(elems(element1(ieli1,1),3)));
            sigma_xy(i)=sigma_xy(i)+0.5*(sigma_xy(elems(element1(ieli1,1),2))+sigma_xy(elems(element1(ieli1,1),3)));            
            elseif element1(ieli1,2)==5
            sigma_xx(i)=sigma_xx(i)+0.5*(sigma_xx(elems(element1(ieli1,1),1))+sigma_xx(elems(element1(ieli1,1),3)));
            sigma_yy(i)=sigma_yy(i)+0.5*(sigma_yy(elems(element1(ieli1,1),1))+sigma_yy(elems(element1(ieli1,1),3)));
            sigma_xy(i)=sigma_xy(i)+0.5*(sigma_xy(elems(element1(ieli1,1),1))+sigma_xy(elems(element1(ieli1,1),3)));
            else
            sigma_xx(i)=sigma_xx(i)+0.5*(sigma_xx(elems(element1(ieli1,1),2))+sigma_xx(elems(element1(ieli1,1),1)));
            sigma_yy(i)=sigma_yy(i)+0.5*(sigma_yy(elems(element1(ieli1,1),2))+sigma_yy(elems(element1(ieli1,1),1)));
            sigma_xy(i)=sigma_xy(i)+0.5*(sigma_xy(elems(element1(ieli1,1),2))+sigma_xy(elems(element1(ieli1,1),1)));                
            end
        end
 sigma_xx(i)=sigma_xx(i)/neli1;
        sigma_yy(i)=sigma_yy(i)/neli1;
        sigma_xy(i)=sigma_xy(i)/neli1;
        sigma(:,:,i)=[sigma_xx(i) sigma_xy(i); sigma_xy(i) sigma_yy(i)];
        [V,D]=eig(sigma(:,:,i));
        sigma_1(i)=-D(1,1);
        sigma_2(i)=-D(2,2);
        sigmavec_1(i,:)=V(:,1)';
        sigmavec_2(i,:)=V(:,2)';       
    else
        for ieli=1:neli
            sigma_xx(i)=sigma_xx(i)+Sigma_xx(element(ieli,1),element(ieli,2));
            sigma_yy(i)=sigma_yy(i)+Sigma_yy(element(ieli,1),element(ieli,2));
            sigma_xy(i)=sigma_xy(i)+Sigma_xy(element(ieli,1),element(ieli,2));


        end
        sigma_xx(i)=sigma_xx(i)/neli;
        sigma_yy(i)=sigma_yy(i)/neli;
        sigma_xy(i)=sigma_xy(i)/neli;
        sigma(:,:,i)=[sigma_xx(i) sigma_xy(i); sigma_xy(i) sigma_yy(i)];
        [V,D]=eig(sigma(:,:,i));
        sigma_1(i)=-D(1,1);
        sigma_2(i)=-D(2,2);
        sigmavec_1(i,:)=-V(:,1)';
        sigmavec_2(i,:)=-V(:,2)';
    cornernodes(i)=i;
    end
    
    
end
toc

cornernodes=nonzeros(cornernodes);
end