%Effect of Roughness on Pseudo-Static Bearing Capacity of Shallow Foundations near Slopes% 
%using the Lower Bound Finite Element Method%

%finite element mesh%elements%nodes%


clc
clear
close all

E=10;
B=2*E;   
a=1*B;

beta=45;
H=80;




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
 % edame mosalas ta H delkhah
 for i=(H/E)+2:(H/E)*2+1
xmo(i,:)=xmo((H/E)+1,:);
ymo(i,:)=(i-1)*E;
xmo2(i,:)=xmo2((H/E)+1,:);
ymo2(i,:)=(2*i-1)*E/2;
 end
 
 
 
 %%
 % baraye zire pei tedad bede
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
 % baraye 2*B
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
 % baraye a
 tedada=7;
xa=zeros(2*(H/E)+1,tedada+1);
ya=zeros(2*(H/E)+1,tedada+1);
xa2=zeros(2*(H/E)+1,tedada+1);
ya2=zeros(2*(H/E)+1,tedada+1);
 for i=1:2*(H/E)+1
    for j=1:tedada+1 
 xa(i,j)=(a/tedada)*(j-1);  
 xa2(i,j)=(a/tedada)/2*(2*j-1); %YE RADIF EZAFE
    end
 ya(i,:)=(i-1)*E;
 ya2(i,:)=(2*i-1)*E/2; %YE RADIF EZAFE
 end 
 
 %%
 %toole zamin samte chap H gereftam
 % tedade eleman ro begire
 tedadezamin=10;
xZA=zeros((H/E)+1,tedadezamin+1);
yZA=zeros((H/E)+1,tedadezamin+1);
xZA2=zeros((H/E)+1,tedadezamin+1);
yZA2=zeros((H/E)+1,tedadezamin+1);
 for i=1:(H/E)+1
    for j=1:tedadezamin+1 
 xZA(i,j)=(H/tedadezamin)*(j-1);  
 xZA2(i,j)=(H/tedadezamin)/2*(2*j-1);    %YE RADIF EZAFE
    end
 yZA(i,:)=(i-1)*E+H;
 yZA2(i,:)=(2*i-1)*E/2+H;   %YE RADIF EZAFE
 end 
 
 %% tabdil
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
plot([Xa(i,1),Xa(i,end)],[[Ya(i,1),Ya(i,end)]],'r-')
plot([Xpei(i,1),Xpei(i,end)],[[Ypei(i,1),Ypei(i,end)]],'r-')
plot([-xez(i,1),-xez(i,end)],[[-yez(i,1),-yez(i,end)]],'r-')
hold on
end
axis([-300 40 -200 40])
for i=1:(H/E)+1
plot([Xmo(i,1),Xmo(i,i)],[[Ymo(i,2),Ymo(i,2)]],'r-')
plot([XZA(i,1),XZA(i,end)],[[YZA(i,2),YZA(i,end)]],'r-')
end
for i=(H/E)+1:2*(H/E)+1
plot([Xmo(i,1),Xmo(i,end-1)],[[Ymo(i,2),Ymo(i,end-1)]],'r-')
end
plot([Xmo(1,1),Xmo((H/E)+1,end-1)],[[Ymo(1,2),Ymo((H/E)+1,end-1)]],'r-')

rectangle('Position',[-3*B,0,B,0.2*B],'FaceColor',[0 .5 .5],'EdgeColor','b','LineWidth',2)

for i=1:tedada+1
plot([Xa(1,i),Xa(end,i)],[[Ya(1,i),Ya(end,i)]],'r-')

end
for i=1:tedadpei+1
plot([Xpei(1,i),Xpei(end,i)],[[Ypei(1,i),Ypei(end,i)]],'r-')
end
for i=1:tedadez+1
plot([-xez(1,i),-xez(end,i)],[[-yez(1,i),-yez(end,i)]],'r-')
end
for i=1:tedadezamin+1
plot([XZA(1,i),XZA(end,i)],[[YZA(1,i),YZA(end,i)]],'r-')
end
for i=1:(H/E)
plot([Xmo(i,i),Xmo(end,i)],[[Ymo(i,i),Ymo(end,i)]],'r-')
end


%% rasm zarbdari
for i=1:2*(H/E)
for j=1:tedada
plot([Xa(i,j),Xa(i+1,j+1)],[[Ya(i,j),Ya(i+1,j+1)]],'r-')
plot([Xa(i,j+1),Xa(i+1,j)],[[Ya(i,j+1),Ya(i+1,j)]],'r-')
end
end

for i=1:2*(H/E)
for j=1:tedadez
plot([-xez(i,j),-xez(i+1,j+1)],[[-yez(i,j),-yez(i+1,j+1)]],'r-')
plot([-xez(i,j+1),-xez(i+1,j)],[[-yez(i,j+1),-yez(i+1,j)]],'r-')
end
end

for i=1:2*(H/E)
for j=1:tedadpei
plot([Xpei(i,j),Xpei(i+1,j+1)],[[Ypei(i,j),Ypei(i+1,j+1)]],'r-')
plot([Xpei(i,j+1),Xpei(i+1,j)],[[Ypei(i,j+1),Ypei(i+1,j)]],'r-')
end
end

for i=1:(H/E)
for j=1:tedadezamin
plot([XZA(i,j),XZA(i+1,j+1)],[[YZA(i,j),YZA(i+1,j+1)]],'r-')
plot([XZA(i,j+1),XZA(i+1,j)],[[YZA(i,j+1),YZA(i+1,j)]],'r-')
end
end

for i=(H/E)+1:2*(H/E)
for j=1:(H/E)
plot([Xmo(i,j),Xmo(i+1,j+1)],[[Ymo(i,j),Ymo(i+1,j+1)]],'r-')
plot([Xmo(i,j+1),Xmo(i+1,j)],[[Ymo(i,j+1),Ymo(i+1,j)]],'r-')
end
end

for i=1:(H/E)
for j=1:i-1
plot([Xmo(i,j),Xmo(i+1,j+1)],[[Ymo(i,j),Ymo(i+1,j+1)]],'r-')
plot([Xmo(i,j+1),Xmo(i+1,j)],[[Ymo(i,j+1),Ymo(i+1,j)]],'r-')
end
end
%%

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

plot(x,y,'*')
for i=1:2*(H/E)+1
%y(2*i-1,tedada+tedadpei+tedadez+1)=-E*i;
end
% Tashkhis tedade setr va sotoone har laye
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
%l=1;

%for i=2:n-1
%if (xy(i,2)-xy((i-1),2))==0
%else
         %   set(l)=i;
         %   l=l+1;
%end
%end





for i=1:n-1
 %plot(xy(i,1),xy(i,2),'o')
 STR={i};
text(xy(i,1),xy(i,2),STR)
 hold on
end
for i=1:f-1
 %plot(xy(i,1),xy(i,2),'o')
 STR={i};
text(xyy(i,1),xyy(i,2),STR)
 hold on
end
sotoon=[0,sotoon];
sotoonx=[0,sotoonx];
mm=0;
for i=1:2*(H/E)
    for j=1:sotoon(i+1)-1
        mm=mm+1;
   nn(mm,:)=[sum(sotoon(1,1:i))+j,sum(sotoon(1,1:i))+j+1,sum(sotoonx(1,1:i))+j];
    end
end
mm=0;
for i=2:2*(H/E)+1
    for j=1:sotoon(i+1)-1
        mm=mm+1;
   nn2(mm,:)=[sum(sotoon(1,1:i))+j,sum(sotoon(1,1:i))+j+1,sum(sotoonx(1,1:i-1))+j];
    end
end
mm=0;
for i=1:2*(H/E)
    for j=1:sotoon(i+1)-1
        mm=mm+1;
   nnr(mm,:)=[sum(sotoon(1,1:i))+j,sum(sotoonx(1,1:i))+j,sum(sotoon(1,1:i+1))+j];
    end
end
mm=0;
for i=1:2*(H/E)
    for j=1:sotoon(i+1)-1
        mm=mm+1;
   nnl(mm,:)=[sum(sotoon(1,1:i))+j+1,sum(sotoon(1,1:i+1))+j+1,sum(sotoonx(1,1:i))+j];
    end
end
NUM=length(nnl)+length(nnr)+length(nn2)+length(nn)+H/E;
STR={NUM};
text(1.5,1,STR)
text(0.5*B,1,'ELEMENTS')
for i=2:H/E+1
nez(i-1,:)=[sum(sotoon(1:i)),sum(sotoon(1:i+1)),sum(sotoon(1:i+1))-1];
end
nez(end,2)=nez(end,2)-tedadezamin;
nez(end,3)=nez(end,3)-tedadezamin;
nm=zeros(length(nn)+1+length(nn2)+1,3);
nm(1:length(nn),:)=nn;
m=0;
for i=length(nn)+1:length(nn)+length(nn2)
    m=m+1;
nm(i,:)=nn2(m,:);
end
m=0;
for i=length(nn)+length(nn2)+1:length(nn)+length(nn2)+length(nnr)
    m=m+1;
nm(i,:)=nnr(m,:);
end
m=0;
for i=length(nn)+length(nn2)+1++length(nnr):length(nn)+length(nn2)+2*length(nnr)
    m=m+1;
nm(i,:)=nnl(m,:);
end
nm(length(nn)+length(nn2)+2*length(nnr)+1:length(nn)+length(nn2)+2*length(nnr)+4,:)=nez;