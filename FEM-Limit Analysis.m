%Effect of Roughness on Pseudo-Static Bearing Capacity of Shallow Foundations near Slopes% 
%using the Lower Bound Finite Element Method%

clc
clear
close all

%soil and Foundation properties%
E=10;     
B=2*E;    
a=1*B;    
beta=45;                                                                
H=80;    
phi=20;   
c=1;      
p=12;     
ts=0;     
g=20;     
kh=0.1;    
kv=kh/2;  
qbar=0;    
delta=phi;  
teta=atan(kh/1-kv); 

%construct mesh using triangle elements
for i=1:(H/E)+1
    lenlaye(i)=(E*(i-1))/tand(beta)+a+B+2*B;
    lenmosalas(i)=(E*(i-1))/tand(beta);
    laye(i)=lenlaye(i)/E;
end
bmm=0;
minl=min(lenmosalas(2:end));
for i=minl:-1:2
    j=0;
    for m=2:(H/E)+1
        if rem(lenmosalas(m),i)==0
            j=j+1;
        end
    end
    if j==length(lenmosalas)-1
        bmm=i;
        break
    end
end

xmo=zeros((H/E)+1,(H/E)+1);
ymo=zeros((H/E)+1,(H/E)+1);
xmo2=zeros((H/E)+1,(H/E)+1);
ymo2=zeros((H/E)+1,(H/E)+1);
LAYE=floor(laye);
for i=2:(H/E)+1
     for m=1:(lenmosalas(i)/bmm)
xmo(i,m)=(lenmosalas(i)-bmm*((lenmosalas(i)/bmm)-m));
xmo2(i,m)=(bmm/2)*(2*m-1);
ymo(i,m)=(i-1)*E;
ymo2(i,m)=(2*i-1)*E/2;
     end
end
 % under slope region
 for i=(H/E)+2:(H/E)*2+1
xmo(i,:)=xmo((H/E)+1,:);
ymo(i,:)=(i-1)*E;
xmo2(i,:)=xmo2((H/E)+1,:);
ymo2(i,:)=(2*i-1)*E/2;
 end
 
 
 
 %%
 % under foundation region
tedadpei=30;
xpei=zeros(2*(H/E)+1,tedadpei+1);
ypei=zeros(2*(H/E)+1,tedadpei+1);
xpei2=zeros(2*(H/E)+1,tedadpei+1);
ypei2=zeros(2*(H/E)+1,tedadpei+1);
 for i=1:2*(H/E)+1
    for j=1:tedadpei+1
 xpei(i,j)=(B/tedadpei)*(j-1);  
 xpei2(i,j)=(B/tedadpei)/2*(2*j-1);
    end
 ypei(i,:)=(i-1)*E;
 ypei2(i,:)=(2*i-1)*E/2;
 end
 %%
% right hand of foundation, assume that 2*B
tedadez=20;
xez=zeros(2*(H/E)+1,tedadez+1);
yez=zeros(2*(H/E)+1,tedadez+1);
xez2=zeros(2*(H/E)+1,tedadez+1);
yez2=zeros(2*(H/E)+1,tedadez+1);
 for i=1:2*(H/E)+1
    for j=1:tedadez+1
 xez(i,j)=(2*B/tedadez)*(j-1);  
 xez2(i,j)=(2*B/tedadez)/2*(2*j-1);
    end
 yez(i,:)=(i-1)*E;
 yez2(i,:)=(2*i-1)*E/2;
 end 
 
 %%
 % distance between slope and foundation region
 tedada=7;
xa=zeros(2*(H/E)+1,tedada+1);
ya=zeros(2*(H/E)+1,tedada+1);
xa2=zeros(2*(H/E)+1,tedada+1);
ya2=zeros(2*(H/E)+1,tedada+1);
 for i=1:2*(H/E)+1
    for j=1:tedada+1 
 xa(i,j)=(a/tedada)*(j-1);  
 xa2(i,j)=(a/tedada)/2*(2*j-1); %one row exceed
    end
 ya(i,:)=(i-1)*E;
 ya2(i,:)=(2*i-1)*E/2; %one row exceed
 end 
 
 %%
 %left hand of slope, assume that H
 tedadezamin=10;
xZA=zeros((H/E)+1,tedadezamin+1);
yZA=zeros((H/E)+1,tedadezamin+1);
xZA2=zeros((H/E)+1,tedadezamin+1);
yZA2=zeros((H/E)+1,tedadezamin+1);
 for i=1:(H/E)+1
    for j=1:tedadezamin+1 
 xZA(i,j)=(H/tedadezamin)*(j-1);  
 xZA2(i,j)=(H/tedadezamin)/2*(2*j-1);    %one row exceed
    end
 yZA(i,:)=(i-1)*E+H;
 yZA2(i,:)=(2*i-1)*E/2+H;   %one row exceed
 end 
 
 %% 
 %Transmission all regions to the global coordinates
Xmo(1:2*(H/E)+1,1)=-(a+B+2*B);
Xmo(1:2*(H/E)+1,2:(H/E)+2)=-(xmo+a+B+2*B);
Ymo(1:2*(H/E)+1,1)=-ymo(:,1);
Ymo(1:2*(H/E)+1,2:(H/E)+2)=-(ymo);
Xmo2(1:2*(H/E)+1,1)=-(a+B+2*B);
Xmo2(1:2*(H/E)+1,2:(H/E)+2)=-(xmo2+a+B+2*B);
Ymo2(1:2*(H/E)+1,1)=-ymo2(:,1);
Ymo2(1:2*(H/E)+1,2:(H/E)+2)=-(ymo2);
 
Xa=-(xa+B+2*B);
Ya=-(ya);
Xa2=-(xa2+B+2*B);
Ya2=-(ya2);

Xpei=-(xpei+2*B);
Ypei=-(ypei);
Xpei2=-(xpei2+2*B);
Ypei2=-(ypei2);

XZA=-(xZA+3*B+a+H/tand(beta));
YZA=-(yZA);
XZA2=-(xZA2+3*B+a+H/tand(beta));
YZA2=-(yZA2);


%% plot
for i=1:2*(H/E)+1
    hold on
plot([Xa(i,1),Xa(i,end)],[[Ya(i,1),Ya(i,end)]],'k-')
plot([Xpei(i,1),Xpei(i,end)],[[Ypei(i,1),Ypei(i,end)]],'k-')
plot([-xez(i,1),-xez(i,end)],[[-yez(i,1),-yez(i,end)]],'k-')
hold on
end
axis([-300 40 -200 40])
for i=1:(H/E)+1
plot([Xmo(i,1),Xmo(i,i)],[[Ymo(i,2),Ymo(i,2)]],'k-')
plot([XZA(i,1),XZA(i,end)],[[YZA(i,2),YZA(i,end)]],'k-')
end
for i=(H/E)+1:2*(H/E)+1
plot([Xmo(i,1),Xmo(i,end-1)],[[Ymo(i,2),Ymo(i,end-1)]],'k-')
end
plot([Xmo(1,1),Xmo((H/E)+1,end-1)],[[Ymo(1,2),Ymo((H/E)+1,end-1)]],'k-')

rectangle('Position',[-3*B,0,B,0.2*B],'FaceColor',[0 0 0],'EdgeColor','k','LineWidth',2)

for i=1:tedada+1
plot([Xa(1,i),Xa(end,i)],[[Ya(1,i),Ya(end,i)]],'k-')

end
for i=1:tedadpei+1
plot([Xpei(1,i),Xpei(end,i)],[[Ypei(1,i),Ypei(end,i)]],'k-')
end
for i=1:tedadez+1
plot([-xez(1,i),-xez(end,i)],[[-yez(1,i),-yez(end,i)]],'k-')
end
for i=1:tedadezamin+1
plot([XZA(1,i),XZA(end,i)],[[YZA(1,i),YZA(end,i)]],'k-')
end
for i=1:(H/E)
plot([Xmo(i,i),Xmo(end,i)],[[Ymo(i,i),Ymo(end,i)]],'k-')
end


%% 
%crosswise draw
for i=1:2*(H/E)
for j=1:tedada
plot([Xa(i,j),Xa(i+1,j+1)],[[Ya(i,j),Ya(i+1,j+1)]],'k-')
plot([Xa(i,j+1),Xa(i+1,j)],[[Ya(i,j+1),Ya(i+1,j)]],'k-')
end
end

for i=1:2*(H/E)
for j=1:tedadez
plot([-xez(i,j),-xez(i+1,j+1)],[[-yez(i,j),-yez(i+1,j+1)]],'k-')
plot([-xez(i,j+1),-xez(i+1,j)],[[-yez(i,j+1),-yez(i+1,j)]],'k-')
end
end

for i=1:2*(H/E)
for j=1:tedadpei
plot([Xpei(i,j),Xpei(i+1,j+1)],[[Ypei(i,j),Ypei(i+1,j+1)]],'k-')
plot([Xpei(i,j+1),Xpei(i+1,j)],[[Ypei(i,j+1),Ypei(i+1,j)]],'k-')
end
end

for i=1:(H/E)
for j=1:tedadezamin
plot([XZA(i,j),XZA(i+1,j+1)],[[YZA(i,j),YZA(i+1,j+1)]],'k-')
plot([XZA(i,j+1),XZA(i+1,j)],[[YZA(i,j+1),YZA(i+1,j)]],'k-')
end
end

for i=(H/E)+1:2*(H/E)
for j=1:(H/E)
plot([Xmo(i,j),Xmo(i+1,j+1)],[[Ymo(i,j),Ymo(i+1,j+1)]],'k-')
plot([Xmo(i,j+1),Xmo(i+1,j)],[[Ymo(i,j+1),Ymo(i+1,j)]],'k-')
end
end

for i=1:(H/E)
for j=1:i-1
plot([Xmo(i,j),Xmo(i+1,j+1)],[[Ymo(i,j),Ymo(i+1,j+1)]],'k-')
plot([Xmo(i,j+1),Xmo(i+1,j)],[[Ymo(i,j+1),Ymo(i+1,j)]],'k-')
end
end
%%
%Number of triangle elements%
xzamin=zeros(2*(H/E)+1,tedadezamin+1);
xzamin((H/E)+1:2*(H/E)+1,:)=XZA;
xzamin2=zeros(2*(H/E)+1,tedadezamin+1);
xzamin2((H/E)+1:2*(H/E)+1,:)=XZA2;
yzamin=zeros(2*(H/E)+1,tedadezamin+1);
yzamin((H/E)+1:2*(H/E)+1,:)=YZA;
yzamin2=zeros(2*(H/E)+1,tedadezamin+1);
yzamin2((H/E)+1:2*(H/E)+1,:)=YZA2;
for i=1:2*(H/E)+1
x(i,:)=[-xez(i,1:end-1),Xpei(i,1:end-1),Xa(i,1:end-1),Xmo(i,1:end-2),xzamin(i,1:end)];
xx(i,:)=[-xez2(i,1:end-1),Xpei2(i,1:end-1),Xa2(i,1:end-1),Xmo2(i,2:end-1),xzamin2(i,1:end)];
end
for i=1:2*(H/E)+1
y(i,:)=[-yez(i,1:end-1),Ypei(i,1:end-1),Ya(i,1:end-1),Ymo(i,1:end-2),yzamin(i,1:end)];
yy(i,:)=[-yez2(i,1:end-1),Ypei2(i,1:end-1),Ya2(i,1:end-1),Ymo2(i,2:end-1),yzamin2(i,1:end)];
end

%plot(x,y,'*')
for i=1:2*(H/E)+1
%y(2*i-1,tedada+tedadpei+tedadez+1)=-E*i;
end
% local node and element number
for i=1:(H/E)+1
    lenlaye(i)=(E*(i-1))/tand(beta)+a+B+2*B;
    lenmosalas(i)=(E*(i-1))/tand(beta);
    laye(i)=lenlaye(i)/E;
end
for i=1:(H/E)
sotoon(i)=tedadez+tedada+tedadpei+(lenmosalas(i)/bmm)+1;
sotoonx(i)=tedadez+tedada+tedadpei+(lenmosalas(i)/bmm);
end
for i=(H/E)+1:2*(H/E)+1
sotoon(i)=tedadezamin+(lenmosalas((H/E))/bmm)+tedadez+tedada+tedadpei+2;
sotoonx(i)=tedadezamin+(lenmosalas((H/E))/bmm)+tedadez+tedada+tedadpei+1;
end
f=1;
for i=1:2*(H/E)
    for j=1:sotoonx(i)
        xyy(f,1)=xx(i,j);
        xyy(f,2)=yy(i,j);
        f=f+1;
       
    end
end

n=1;
for i=1:2*(H/E)+1
    for j=1:sotoon(i)
        xy(n,1)=x(i,j);
        xy(n,2)=y(i,j);
        n=n+1;     
     
    end
end
l=1;

for i=2:n-1
if (xy(i,2)-xy((i-1),2))==0
else
           set(l)=i;
            l=l+1;
end
end




%STR={n-1};
%text(1,10,STR)
%text(0.5*B,10,'nodes')
%for i=1:n-1%text of node numbers
%plot(xy(i,1),xy(i,2),'o')
%STR={i};%text of node numbers
%text(xy(i,1),xy(i,2),STR)%text of node numbers
%hold on%text of node numbers
%end%text of node numbers
%for i=1:f-1%text of node numbers
 %plot(xy(i,1),xy(i,2),'o')  %hazf she
 %STR={i};%text of node numbers
%text(xyy(i,1),xyy(i,2),STR)%text of node numbers
 %hold on%text of node numbers
%end%text of node numbers
sotoon=[0,sotoon];
sotoonx=[0,sotoonx];
%4 type triangle
%Up triangle
mm=0;
for i=1:2*(H/E)
    for j=1:sotoon(i+1)-1
        mm=mm+1;
   nn(mm,:)=[sum(sotoon(1,1:i))+j,sum(sotoon(1,1:i))+j+1,sum(sotoonx(1,1:i))+j];
    end
end
%Down triangle
mm=0;
for i=2:2*(H/E)+1
    for j=1:sotoon(i+1)-1
        mm=mm+1;
   nn2(mm,:)=[sum(sotoon(1,1:i))+j,sum(sotoon(1,1:i))+j+1,sum(sotoonx(1,1:i-1))+j];
    end
end
%Right triangle
mm=0;
for i=1:2*(H/E)
    for j=1:sotoon(i+1)-1
        mm=mm+1;
   nnr(mm,:)=[sum(sotoon(1,1:i))+j,sum(sotoonx(1,1:i))+j,sum(sotoon(1,1:i+1))+j];
    end
end
%Left triangle
mm=0;
for i=1:2*(H/E)
    for j=1:sotoon(i+1)-1
        mm=mm+1;
   nnl(mm,:)=[sum(sotoon(1,1:i))+j+1,sum(sotoon(1,1:i+1))+j+1,sum(sotoonx(1,1:i))+j];
    end
end
%number of the triangle elements
NUM=length(nnl)+length(nnr)+length(nn2)+length(nn)+H/E;
STR={NUM};
text(1.5,1,STR)
text(0.5*B,1,'ELEMENTS')
%Shape Functions

%up triangles shape functions
au=zeros(3,length(nn));bu=zeros(3,length(nn));cu=zeros(3,length(nn));

for i=1:length(nn)
    au(1,i)=(xy(nn(i,2),1)*xyy(nn(i,3),2))-(xyy(nn(i,3),1)*xy(nn(i,2),2));
    au(2,i)=(xyy(nn(i,3),1)*xy(nn(i,1),2))-(xy(nn(i,1),1)*xyy(nn(i,3),2));
    au(3,i)=(xy(nn(i,1),1)*xy(nn(i,2),2))-(xy(nn(i,2),1)*xy(nn(i,1),2));
    bu(1,i)=xy(nn(i,2),2)-xyy(nn(i,3),2);
    bu(2,i)=xyy(nn(i,3),2)-xy(nn(i,1),2);
    bu(3,i)=xy(nn(i,1),2)-xy(nn(i,2),2);
    cu(1,i)=xyy(nn(i,3),1)-xy(nn(i,2),1);
    cu(2,i)=xy(nn(i,1),1)-xyy(nn(i,3),1);
    cu(3,i)=xy(nn(i,2),1)-xy(nn(i,1),1);
end


%down triangles shape functions
ad=zeros(3,length(nn2));bd=zeros(3,length(nn2));cd=zeros(3,length(nn2));

for i=1:length(nn2)
    ad(1,i)=(xy(nn2(i,2),1)*xyy(nn2(i,3),2))-(xyy(nn2(i,3),1)*xy(nn2(i,2),2));
    ad(2,i)=(xyy(nn2(i,3),1)*xy(nn2(i,1),2))-(xy(nn2(i,1),1)*xyy(nn2(i,3),2));
    ad(3,i)=(xy(nn2(i,1),1)*xy(nn2(i,2),2))-(xy(nn2(i,2),1)*xy(nn2(i,1),2));
    bd(1,i)=xy(nn2(i,2),2)-xyy(nn2(i,3),2);
    bd(2,i)=xyy(nn2(i,3),2)-xy(nn2(i,1),2);
    bd(3,i)=xy(nn2(i,1),2)-xy(nn2(i,2),2);
    cd(1,i)=-1*(xyy(nn2(i,3),1)-xy(nn2(i,2),1));
    cd(2,i)=-1*(xy(nn2(i,1),1)-xyy(nn2(i,3),1));
    cd(3,i)=-1*(xy(nn2(i,2),1)-xy(nn2(i,1),1));
end

%right triangles shape functions
ar=zeros(3,length(nnr));br=zeros(3,length(nnr));cr=zeros(3,length(nnr));

for i=1:length(nnr)
    ar(1,i)=(xyy(nnr(i,2),1)*xy(nnr(i,3),2))-(xy(nnr(i,3),1)*xyy(nnr(i,2),2));
    ar(2,i)=(xy(nnr(i,3),1)*xy(nnr(i,1),2))-(xy(nnr(i,1),1)*xy(nnr(i,3),2));
    ar(3,i)=(xy(nnr(i,1),1)*xyy(nnr(i,2),2))-(xyy(nnr(i,2),1)*xy(nnr(i,1),2));
    br(1,i)=xyy(nnr(i,2),2)-xy(nnr(i,3),2);
    br(2,i)=xy(nnr(i,3),2)-xy(nnr(i,1),2);
    br(3,i)=xy(nnr(i,1),2)-xyy(nnr(i,2),2);
    cr(1,i)=xy(nnr(i,3),1)-xyy(nnr(i,2),1);
    cr(2,i)=xy(nnr(i,1),1)-xy(nnr(i,3),1);
    cr(3,i)=xyy(nnr(i,2),1)-xy(nnr(i,1),1);
end


%left triangles shape functions
al=zeros(3,length(nnl));bl=zeros(3,length(nnl));cl=zeros(3,length(nnl));

for i=1:length(nnl)
    al(1,i)=(xy(nnl(i,2),1)*xyy(nnl(i,3),2))-(xyy(nnl(i,3),1)*xy(nnl(i,2),2));
    al(2,i)=(xyy(nnl(i,3),1)*xy(nnl(i,1),2))-(xy(nnl(i,1),1)*xyy(nnl(i,3),2));
    al(3,i)=(xy(nnl(i,1),1)*xy(nnl(i,2),2))-(xy(nnl(i,2),1)*xy(nnl(i,1),2));
    bl(1,i)=xy(nnl(i,2),2)-xyy(nnl(i,3),2);
    bl(2,i)=xyy(nnl(i,3),2)-xy(nnl(i,1),2);
    bl(3,i)=xy(nnl(i,1),2)-xy(nnl(i,2),2);
    cl(1,i)=xyy(nnl(i,3),1)-xy(nnl(i,2),1);
    cl(2,i)=xy(nnl(i,1),1)-xyy(nnl(i,3),1);
    cl(3,i)=xy(nnl(i,2),1)-xy(nnl(i,1),1);
end


%Assembling triangles
%for i=2:H/E+1
%nez(i-1,:)=[sum(sotoon(1:i)),sum(sotoon(1:i+1)),sum(sotoon(1:i+1))-1];
%end
%nez(end,2)=nez(end,2)-tedadezamin;
%nez(end,3)=nez(end,3)-tedadezamin;
%nm=zeros(length(nn)+1+length(nn2)+1,3);
%nm(1:length(nn),:)=nn;
%m=0;
%for i=length(nn)+1:length(nn)+length(nn2)
%   m=m+1;
%nm(i,:)=nn2(m,:);
%end
%m=0;
%for i=length(nn)+length(nn2)+1:length(nn)+length(nn2)+length(nnr)
 %   m=m+1;
%nm(i,:)=nnr(m,:);
%end
%m=0;
%for i=length(nn)+length(nn2)+1++length(nnr):length(nn)+length(nn2)+2*length(nnr)
 %   m=m+1;
%nm(i,:)=nnl(m,:);
%end
%nm(length(nn)+length(nn2)+2*length(nnr)+1:length(nn)+length(nn2)+2*length(nnr)+4,:)=nez;

%Assembling shape functions
 %aa=zeros(3,NUM);bb=zeros(3,NUM);cc=zeros(3,NUM);
 %for i=1:length(nm)
 %   aa(1,i)=ad(1,i)+ad(1,i)+ar(1,i)+al(1,i);
 %   aa(2,i)=ad(2,i)+ad(2,i)+ar(2,i)+al(2,i);
 %   aa(3,i)=ad(3,i)+ad(3,i)+ar(3,i)+al(3,i);
 %   bb(1,i)=bd(1,i)+bd(1,i)+br(1,i)+bl(1,i);
 %   bb(2,i)=bd(2,i)+bd(2,i)+br(2,i)+bl(2,i);
 %  bb(3,i)=bd(3,i)+bd(3,i)+br(3,i)+bl(3,i);
 %  cc(1,i)=cd(1,i)+cd(1,i)+cr(1,i)+cl(1,i);
 %  cc(2,i)=cd(2,i)+cd(2,i)+cr(2,i)+cl(2,i);
 %   cc(3,i)=cd(3,i)+cd(3,i)+cr(3,i)+cl(3,i);
%end

areau=zeros(1,length(nn));                                                            %Area of up triangle Elements.
for i=1:length(nn)
    areau(1,i)=0.5*abs(bu(1,i)*cu(2,i)-bu(2,i)*cu(1,i));
end

aread=zeros(1,length(nn2));                                                            %Area of down triangle Elements.
for i=1:length(nn2)
    aread(1,i)=0.5*abs(bd(1,i)*cd(2,i)-bd(2,i)*cd(1,i));
end

arear=zeros(1,length(nnr));                                                            %Area of right triangle Elements.
for i=1:length(nnr)
    arear(1,i)=0.5*abs(br(1,i)*cr(2,i)-br(2,i)*cr(1,i));
end

areal=zeros(1,length(nnl));                                                            %Area of left triangle Elements.
for i=1:length(nnl)
    areal(1,i)=0.5*abs(bl(1,i)*cl(2,i)-bl(2,i)*cl(1,i));
end
   
%objective Function and Equality & Inequality Constraints

    f=sparse(9*NUM,1);                                   %Objective Function Coeficient Vector.
    A1=sparse(3*p*NUM,9*NUM);                            %Inequality Matrix.
    A2=sparse(2*NUM+4*(NUM-1),9*NUM);                  %Equality Matrix.
    b1=sparse(3*p*NUM,1);
    b2=sparse(2*NUM+4*(NUM-1),1);
    f(1)=0.5*(sind(ts))^2;f(2)=0.5*(cosd(ts))^2;f(3)=-0.5*sind(2*ts);
    f(4)=0.5*(sind(ts))^2;f(5)=0.5*(cosd(ts))^2;f(6)=-0.5*sind(2*ts);
    
%Yield Condition

A=zeros(1,p);BB=zeros(1,p);C=zeros(1,p);                  %Coefficient of Stresses in a Linearized Yield Criterion.
    for kk=1:p
        A(kk)=cos(2*pi*kk/p)+sind(phi)*cos(pi/p);
        BB(kk)=sind(phi)*cos(pi/p)-cos(2*pi*kk/p);
        C(kk)=2*sin(2*pi*kk/p);
    end
    D=2*c*cosd(phi)*cos(pi/p);
    
    K=[A',BB',C'];
    
    m=zeros(1,3*NUM);                              %Some Numerator For Assembling the Matrix.                                                                 
    m(1)=1;
    for ii=2:3*NUM
        m(ii)=m(ii-1)+3;
    end
    
    for jj=1:3*NUM
        A1((jj-1)*p+1:jj*p,m(jj):m(jj)+2)=K;
    end
    
    for ii=1:3*p*NUM
        b1(ii,1)=D;
    end


%Element Equilibrium for up triangles

m=zeros(1,length(nn));w=zeros(1,length(nn));                      %Some Numerator For Assembling the Matrix.
m(1)=1;w(1)=1;
for ii=2:length(nn)
    m(ii)=m(ii-1)+2;w(ii)=w(ii-1)+9;
end

for ii=1:length(nn)
A2(m(ii):m(ii)+1,w(ii):w(ii)+8)=(1/(2*areau(ii)))*[bu(1,ii) 0 cu(1,ii) bu(2,ii) 0 cu(2,ii) bu(3,ii) 0 cu(3,ii);0 cu(1,ii) bu(1,ii) 0 cu(2,ii) bu(2,ii) 0 cu(3,ii) bu(3,ii)];
end

b22=[-kh*g;(1-kv)*g];
for ii=1:length(nn)
    b2(m(ii):m(ii)+1,1)=b22;
end

%Element Equilibrium for down triangles

m=zeros(1,length(nn2));w=zeros(1,length(nn2));                      %Some Numerator For Assembling the Matrix.
m(1)=1;w(1)=1;
for ii=2:length(nn2)
    m(ii)=m(ii-1)+2;w(ii)=w(ii-1)+9;
end

for ii=1:length(nn2)
A2(m(ii)+2:m(ii)+3,w(ii)+9:w(ii)+17)=(1/(2*aread(ii)))*[bd(1,ii) 0 cd(1,ii) bd(2,ii) 0 cd(2,ii) bd(3,ii) 0 cd(3,ii);0 cd(1,ii) bd(1,ii) 0 cd(2,ii) bd(2,ii) 0 cd(3,ii) bd(3,ii)];
end

b22=[-kh*g;(1-kv)*g];
for ii=1:length(nn2)
    b2(m(ii)+2:m(ii)+3,1)=b22;
end

%Element Equilibrium for right triangles

m=zeros(1,length(nnr));w=zeros(1,length(nnr));                      %Some Numerator For Assembling the Matrix.
m(1)=1;w(1)=1;
for ii=2:length(nnr)
    m(ii)=m(ii-1)+2;w(ii)=w(ii-1)+9;
end

for ii=1:length(nnr)
A2(m(ii)+4:m(ii)+5,w(ii)+18:w(ii)+26)=(1/(2*arear(ii)))*[br(1,ii) 0 cr(1,ii) br(2,ii) 0 cr(2,ii) br(3,ii) 0 cr(3,ii);0 cr(1,ii) br(1,ii) 0 cr(2,ii) br(2,ii) 0 cr(3,ii) br(3,ii)];
end

b22=[-kh*g;(1-kv)*g];
for ii=1:length(nnr)
    b2(m(ii)+4:m(ii)+5,1)=b22;
end

%Element Equilibrium for left triangles

m=zeros(1,length(nnl));w=zeros(1,length(nnl));                      %Some Numerator For Assembling the Matrix.
m(1)=1;w(1)=1;
for ii=2:length(nnl)
    m(ii)=m(ii-1)+2;w(ii)=w(ii-1)+9;
end

for ii=1:length(nnl)
A2(m(ii)+6:m(ii)+7,w(ii)+27:w(ii)+35)=(1/(2*areal(ii)))*[bl(1,ii) 0 cl(1,ii) bl(2,ii) 0 cl(2,ii) bl(3,ii) 0 cl(3,ii);0 cl(1,ii) bl(1,ii) 0 cl(2,ii) bl(2,ii) 0 cl(3,ii) bl(3,ii)];
end

b22=[-kh*g;(1-kv)*g];
for ii=1:length(nnl)
    b2(m(ii)+6:m(ii)+7,1)=b22;
end


%disconytinuity lines

for i=1:length(nnr)
    annrf(i)=((xy(nnr(i,1),2)-xyy(nnr(i,2),2))/(xy(nnr(i,1),1)-xyy(nnr(i,2),1)));    %discontinuity line type1
    annrb(i)=((xyy(nnr(i,2),2)-xy(nnr(i,3),2))/(xyy(nnr(i,2),1)-xy(nnr(i,3),1)));    %discontinuity line type2
end

for i=1:length(nnl)
    annlf(i)=((xy(nnl(i,1),2)-xyy(nnl(i,3),2))/(xy(nnl(i,1),1)-xyy(nnl(i,3),1)));     %disconyinuity line type3
    annlb(i)=((xy(nnl(i,2),2)-xyy(nnl(i,3),2))/(xy(nnl(i,2),1)-xyy(nnl(i,3),1)));     %disconyinuity line type4
end

%Discontinuity Equilibrium
%discontinuity line type 1
tdrf=zeros(1,length(nnr)-1);                             %Counter Clockwise Inclination of Each Inner Mesh Sides.

for ii=1:length(nnr)-1
   tdrf(1,ii)=atan(annrf(i));
end

m=zeros(1,length(nnr)-1);w=zeros(1,length(nnr)-1);               %Some Numerator For Assembling the Matrix.
m(1,1)=1;w(1,1)=1;
for ii=2:length(nnr)-1
    m(1,ii)=m(1,ii-1)+4;
    w(1,ii)=w(1,ii-1)+12;
end
for ii=1:length(nnr)-1
    A2((2*NUM+m(ii)):(2*NUM+m(ii)+3),w(ii):(w(ii)+11))=[sin(tdrf(ii))^2 cos(tdrf(ii))^2 -sin(2*tdrf(ii)) -sin(tdrf(ii))^2 -cos(tdrf(ii))^2 sin(2*tdrf(ii)) 0 0 0 0 0 0;-0.5*sin(2*tdrf(ii)) 0.5*sin(2*tdrf(ii)) cos(2*tdrf(ii)) 0.5*sin(2*tdrf(ii)) -0.5*sin(2*tdrf(ii)) -cos(2*tdrf(ii)) 0 0 0 0 0 0;0 0 0 0 0 0 sin(tdrf(ii))^2 cos(tdrf(ii))^2 -sin(2*tdrf(ii)) -sin(tdrf(ii))^2 -cos(tdrf(ii))^2 sin(2*tdrf(ii));0 0 0 0 0 0 -0.5*sin(2*tdrf(ii)) 0.5*sin(2*tdrf(ii)) cos(2*tdrf(ii)) 0.5*sin(2*tdrf(ii)) -0.5*sin(2*tdrf(ii)) -cos(2*tdrf(ii))];
end

%disconyinuity line type 2 
tdrb=zeros(1,length(nnr)-1);                                %Counter Clockwise Inclination of Each Inner Mesh Sides.

for ii=1:length(nnr)-1
   tdrb(1,ii)=atan(annrb(i));
end
m=zeros(1,length(nnr)-1);w=zeros(1,length(nnr)-1);               %Some Numerator For Assembling the Matrix.
m(1,1)=1;w(1,1)=1;
for ii=2:length(nnr)-1
    m(1,ii)=m(1,ii-1)+4;
    w(1,ii)=w(1,ii-1)+12;
end
for ii=1:length(nnr)-1
   A2((2*NUM+m(ii))+4:(2*NUM+m(ii)+7),w(ii)+12:(w(ii)+23))=[sin(tdrb(ii))^2 cos(tdrb(ii))^2 -sin(2*tdrb(ii)) -sin(tdrb(ii))^2 -cos(tdrb(ii))^2 sin(2*tdrb(ii)) 0 0 0 0 0 0;-0.5*sin(2*tdrb(ii)) 0.5*sin(2*tdrb(ii)) cos(2*tdrb(ii)) 0.5*sin(2*tdrb(ii)) -0.5*sin(2*tdrb(ii)) -cos(2*tdrb(ii)) 0 0 0 0 0 0;0 0 0 0 0 0 sin(tdrb(ii))^2 cos(tdrb(ii))^2 -sin(2*tdrb(ii)) -sin(tdrb(ii))^2 -cos(tdrb(ii))^2 sin(2*tdrb(ii));0 0 0 0 0 0 -0.5*sin(2*tdrb(ii)) 0.5*sin(2*tdrb(ii)) cos(2*tdrb(ii)) 0.5*sin(2*tdrb(ii)) -0.5*sin(2*tdrb(ii)) -cos(2*tdrb(ii))]; 
end

%disconyinuity line type 3 
tdlf=zeros(1,length(nnl)-1);                             %Counter Clockwise Inclination of Each Inner Mesh Sides.

for ii=1:length(nnl)-1
   tdlf(1,ii)=atan(annlf(i));
end
m=zeros(1,length(nnl)-1);w=zeros(1,length(nnl)-1);               %Some Numerator For Assembling the Matrix.
m(1,1)=1;w(1,1)=1;
for ii=2:length(nnl)-1
    m(1,ii)=m(1,ii-1)+4;
    w(1,ii)=w(1,ii-1)+12;
end
for ii=1:length(nnl)-1
     A2((2*NUM+m(ii)+4):(2*NUM+m(ii)+7),w(ii)+24:(w(ii)+35))=[sin(tdlf(ii))^2 cos(tdlf(ii))^2 -sin(2*tdlf(ii)) -sin(tdlf(ii))^2 -cos(tdlf(ii))^2 sin(2*tdlf(ii)) 0 0 0 0 0 0;-0.5*sin(2*tdlf(ii)) 0.5*sin(2*tdlf(ii)) cos(2*tdlf(ii)) 0.5*sin(2*tdlf(ii)) -0.5*sin(2*tdlf(ii)) -cos(2*tdlf(ii)) 0 0 0 0 0 0;0 0 0 0 0 0 sin(tdlf(ii))^2 cos(tdlf(ii))^2 -sin(2*tdlf(ii)) -sin(tdlf(ii))^2 -cos(tdlf(ii))^2 sin(2*tdlf(ii));0 0 0 0 0 0 -0.5*sin(2*tdlf(ii)) 0.5*sin(2*tdlf(ii)) cos(2*tdlf(ii)) 0.5*sin(2*tdlf(ii)) -0.5*sin(2*tdlf(ii)) -cos(2*tdlf(ii))];
end

%disconyinuity line type 4 
tdlb=zeros(1,length(nnl)-1);                             %Counter Clockwise Inclination of Each Inner Mesh Sides.

for ii=1:length(nnl)-1
   tdlb(1,ii)=atan(annlb(i));
end
m=zeros(1,length(nnl)-1);w=zeros(1,length(nnl)-1);               %Some Numerator For Assembling the Matrix.
m(1,1)=1;w(1,1)=1;
for ii=2:length(nnl)-1
    m(1,ii)=m(1,ii-1)+4;
    w(1,ii)=w(1,ii-1)+12;
end
for ii=1:length(nnl)-1
    A2((2*NUM+m(ii)+8):(2*NUM+m(ii)+11),w(ii)+36:(w(ii)+47))=[sin(tdlb(ii))^2 cos(tdlb(ii))^2 -sin(2*tdlb(ii)) -sin(tdlb(ii))^2 -cos(tdlb(ii))^2 sin(2*tdlb(ii)) 0 0 0 0 0 0;-0.5*sin(2*tdlb(ii)) 0.5*sin(2*tdlb(ii)) cos(2*tdlb(ii)) 0.5*sin(2*tdlb(ii)) -0.5*sin(2*tdlb(ii)) -cos(2*tdlb(ii)) 0 0 0 0 0 0;0 0 0 0 0 0 sin(tdlb(ii))^2 cos(tdlb(ii))^2 -sin(2*tdlb(ii)) -sin(tdlb(ii))^2 -cos(tdlb(ii))^2 sin(2*tdlb(ii));0 0 0 0 0 0 -0.5*sin(2*tdlb(ii)) 0.5*sin(2*tdlb(ii)) cos(2*tdlb(ii)) 0.5*sin(2*tdlb(ii)) -0.5*sin(2*tdlb(ii)) -cos(2*tdlb(ii))];
end

%discontinuity line type 5 %horizontol discontinuity
tduz=zeros(1,length(nn)-1);                             %Counter Clockwise Inclination of Each Inner Mesh Sides.

for ii=1:length(nn)-1
   tduz(1,ii)=0;
end

m=zeros(1,length(nn)-1);w=zeros(1,length(nn)-1);               %Some Numerator For Assembling the Matrix.
m(1,1)=1;w(1,1)=1;
for ii=2:length(nn)-1
    m(1,ii)=m(1,ii-1)+4;
    w(1,ii)=w(1,ii-1)+12;
end
for ii=1:length(nn)-1
     A2((2*NUM+m(ii)+12):(2*NUM+m(ii)+15),w(ii)+48:(w(ii)+59))=[sin(tduz(ii))^2 cos(tduz(ii))^2 -sin(2*tduz(ii)) -sin(tduz(ii))^2 -cos(tduz(ii))^2 sin(2*tduz(ii)) 0 0 0 0 0 0;-0.5*sin(2*tduz(ii)) 0.5*sin(2*tduz(ii)) cos(2*tduz(ii)) 0.5*sin(2*tduz(ii)) -0.5*sin(2*tduz(ii)) -cos(2*tduz(ii)) 0 0 0 0 0 0;0 0 0 0 0 0 sin(tduz(ii))^2 cos(tduz(ii))^2 -sin(2*tduz(ii)) -sin(tduz(ii))^2 -cos(tduz(ii))^2 sin(2*tduz(ii));0 0 0 0 0 0 -0.5*sin(2*tduz(ii)) 0.5*sin(2*tduz(ii)) cos(2*tduz(ii)) 0.5*sin(2*tduz(ii)) -0.5*sin(2*tduz(ii)) -cos(2*tduz(ii))];
end

%discontinuity line type 6 %vertical discontinuity
tdlv=zeros(1,length(nnl)-1);                             %Counter Clockwise Inclination of Each Inner Mesh Sides.

for ii=1:length(nnl)-1
   tdlv(1,ii)=90;
end

m=zeros(1,length(nnl)-1);w=zeros(1,length(nnl)-1);               %Some Numerator For Assembling the Matrix.
m(1,1)=1;w(1,1)=1;
for ii=2:length(nnl)-1
    m(1,ii)=m(1,ii-1)+4;
    w(1,ii)=w(1,ii-1)+12;
end
for ii=1:length(nnl)-1
    A2((2*NUM+m(ii)+16):(2*NUM+m(ii)+19),w(ii)+60:(w(ii)+71))=[sin(tdlv(ii))^2 cos(tdlv(ii))^2 -sin(2*tdlv(ii)) -sin(tdlv(ii))^2 -cos(tdlv(ii))^2 sin(2*tdlv(ii)) 0 0 0 0 0 0;-0.5*sin(2*tdlv(ii)) 0.5*sin(2*tdlv(ii)) cos(2*tdlv(ii)) 0.5*sin(2*tdlv(ii)) -0.5*sin(2*tdlv(ii)) -cos(2*tdlv(ii)) 0 0 0 0 0 0;0 0 0 0 0 0 sin(tdlv(ii))^2 cos(tdlv(ii))^2 -sin(2*tdlv(ii)) -sin(tdlv(ii))^2 -cos(tdlv(ii))^2 sin(2*tdlv(ii));0 0 0 0 0 0 -0.5*sin(2*tdlv(ii)) 0.5*sin(2*tdlv(ii)) cos(2*tdlv(ii)) 0.5*sin(2*tdlv(ii)) -0.5*sin(2*tdlv(ii)) -cos(2*tdlv(ii))];
end
%Boundary Condition

m=zeros(1,length(nn)-1);w=zeros(1,length(nn)-1);               %Some Numerator For Assembling the Matrix.
m(1,1)=1;w(1,1)=1;
for ii=2:length(nn)-1
    m(1,ii)=m(1,ii-1)+4;
    w(1,ii)=w(1,ii-1)+12;
end
for ii=1:10
     A2((2*NUM+m(ii)+12):(2*NUM+m(ii)+15),w(ii)+48:(w(ii)+59))=1;  
end
 
for ii=11:20
     A2((2*NUM+m(ii)+12):(2*NUM+m(ii)+15),w(ii)+48:(w(ii)+59))=1;  
end

for ii=21:30
     A2((2*NUM+m(ii)+12):(2*NUM+m(ii)+15),w(ii)+48:(w(ii)+59))=1;   
end

for ii=165:175
     A2((2*NUM+m(ii)+12):(2*NUM+m(ii)+15),w(ii)+48:(w(ii)+59))=1;  
end

for ii=11:20
     b2((2*NUM+m(ii)+12):(2*NUM+m(ii)+15),1)=-qbar;  
end

%Linear Programming

ss=linprog(-f,A1,b1,A2,b2);

%The Answers
format bank
qu=abs(f'*ss)                                                                         %The Ultimate Bearing Capacity (kPa)
Ng=abs(qu)/(0.5*1*g)


 
