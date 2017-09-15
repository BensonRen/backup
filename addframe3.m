data1 = dlmread('1.dat');
data2 = dlmread('2.dat');

[h1, l1] = size(data1);

m.d=2;                      %dimension
m.nframe1=max(data1(:,l1-1))+1;
m.n1=max(data1(:,l1))+1;            %total number of particle


[h2, l2] = size(data2);
m.nframe2=max(data2(:,l2-1))+1;     %total number of frame of the last particle
m.n2=max(data2(h2,l2))+1;            %total number of particle
m.T = ones(m.nframe1+m.nframe2, 1) * 1;
m.nframe=m.nframe1+m.nframe2;
%m.n=max(m.n1,m.n2);
m.n=m.n1+m.n2;


%sort
%m.radius1 = zeros(m.nframe1,m.n1);
m.interp=ones(m.nframe,m.n);
m.radius=zeros(m.nframe,m.n);
m.r=zeros(m.nframe,m.n,m.d);x=0;
for i=1:h1;
    m.radius(data1(i,l1-1)+1,data1(i,l1)+1) = data1(i,4);
    m.r(data1(i,l1-1)+1,data1(i,l1)+1,1)=data1(i,1);
    m.r(data1(i,l1-1)+1,data1(i,l1)+1,2)=data1(i,2);
    m.interp(data1(i,l1-1)+1,data1(i,l1)+1)=0;
end;


%sort
m.radius2 = zeros(m.nframe2,m.n2);
m.r2=zeros(m.nframe2,m.n2,m.d);    %position of particles in different frames
m.interp2=ones(m.nframe2,m.n2);
for i=1:h2;
    m.radius2(data2(i,l2-1)+1,data2(i,l2)+1) = data2(i,4);
    m.r2(data2(i,l2-1)+1,data2(i,l2)+1,1)=data2(i,1);
    m.r2(data2(i,l2-1)+1,data2(i,l2)+1,2)=data2(i,2);
    %m.r(data1(i,l1-1)+1+h1,data1(i,l1)+1,1)=data2(i,1);
    %m.r(data1(i,l1-1)+1+h1,data1(i,l1)+1,2)=data2(i,2);
    m.interp2(data2(i,l2-1)+1,data2(i,l2)+1)=0;
end;
m.r(:,:,3)=m.r(:,:,1)*0;
m.r2(:,:,3)=m.r2(:,:,1)*0;
rmax=squeeze(max(max(m.r2,[],1),[],2));
rmin=squeeze(min(min(m.r2,[],1),[],2));
L=(rmax-rmin);
diameter=(L(1)*L(2)/m.n2).^0.5;

r1last=squeeze(m.r(m.nframe1,:,:));
r2first=squeeze(m.r2(1,:,:));
drmin=diameter;y=0;u=0;
%drlf=zeros(m.n2);
for i=1:1:m.n
    drmin=diameter;
    if m.interp(m.nframe1,i)==0
        b=0;
        for a=1:1:m.n2
            x=any(a==y);
        if m.interp2(1,a)==0&&x==0
        
            drlf=((r2first(a,1)-r1last(i,1))^2+(r2first(a,2)-r1last(i,2))^2)^0.5;
            if drlf<drmin
               drmin=drlf;
               b=a;
               u=u+1;
            end
              
        end
        end
        
            if b~=0
                if any(y==b)==0
                    
                   for k=1:1:m.nframe2
                       if m.interp2(k,b)==0
                       m.r(m.nframe1+k,i,:)=m.r2(k,b,:); 
                       %m.r(m.nframe1+k,i,2)=m.r2(k,b,2); 
                       %m.interp(m.nframe1+k,i)=0;
                       m.interp(m.nframe1+k,i)=0;
                       m.radius(m.nframe1+k,i) = data2(b,4);
                       else
                        m.interp(m.nframe1+k,i)=1;   
                       
                       end
                   end
                end
                 y=[y;b];
            else
                
                 
            end
     end
            
            
end
    

y(1)=[];

          
u
y
k=0
for i=1:1:m.n2
    if m.interp2(1,i)~=0&&any(i~=y)
        k=[k;i];
    end
end

[p Q]=size(k);

t=max(y)+1

 for n=2:1:(p+1)
     
    if x==1
        for k=1:1:m.nframe2
         if m.interp2(k,p(n,1))==0
            m.r(m.nframe1+k,t,:)=m.r2(k,p(n,1),:); 
            m.interp(m.nframe1+k,t)=0;
            m.radius(m.nframe1+k,t) = data1(p(n,1),4);
         
         end
        end
        t=t+1;
    end
 end
 t
 c=zeros(h1+h2,6);
 
 z=1; 
 
  for j=1:1:m.n
                     
         for i=1:m.nframe
             
             y=m.interp(i,j);
             if y==0
             c(z,6)=j;
             c(z,5)=i;
             c(z,1)=m.r(i,j,1);
             c(z,2)=m.r(i,j,2);
             c(z,4)=m.radius(i,j);
             z=z+1;
             end
         end
       
  end
  if z~=1
      z=z-1;
  end
  
  c=c(1:(z-1),:);
  c(:,5)=c(:,5)-1;c(:,6)=c(:,6)-1;
  [h,l]=size(c)
  
  z
  
 
   z=0;a=c(1,6);

for j=1:1:h
    if a~=c(j,6)
        z=z+1;a=c(j,6);
  
    end
  c(j,6)=z;
  
end

max(c(:,6)) 
 [h,l]=size(c)
 'ok'

 %for i=1:1:h
  %   if c(i,5)<0
   %  c(i,:)=[];
 %    end
 %end
     
 
 
dlmwrite('c@ Yip.dat', c,'\t');
%  save c.dat c -ascii;
%  exit;