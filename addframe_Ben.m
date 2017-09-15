% %%This addframe3 script add up two .dat file to one c.dat file so that the
% %%matlab would do a larger set of data instead of one set at a time
data1 = dlmread('1.dat');   %read data1 from 1.dat
m1.size=size(data1,1);          %size of data
m1.nf=max(data1(:,5))+1;        %total number of frames
m1.np=max(data1(:,6))+1;        %total number of particle
m1.interp=ones(m1.nf,m1.np);    %pre-allocation for interp (1 for obtained by other methods and 0 for real)
m1.radius=zeros(m1.nf,m1.np);   %pre-allocation of radius
m1.r=zeros(m1.nf,m1.np,2);     %pre-allocation of position
m1.mass=zeros(m1.nf,m1.np);     %pre-allocation of mass
%assign all the values into m1
for i=1:m1.size
    m1.radius(data1(i,5)+1,data1(i,6)+1) = data1(i,4);
    m1.r(data1(i,5)+1,data1(i,6)+1,1)=data1(i,1);
    m1.r(data1(i,5)+1,data1(i,6)+1,2)=data1(i,2);
    m1.interp(data1(i,5)+1,data1(i,6)+1)=0;
    m1.mass(data1(i,5)+1,data1(i,6)+1) = data1(i,3);
end
data1=[];                       %release the RAM of data1 occupied
'read data1 ok'

data2 = dlmread('2.dat');       %read data2 from 2.dat
m2.size=size(data2,1);          %size of data
m2.nf=max(data2(:,5))+1;        %total number of frames
m2.np=max(data2(:,6))+1;        %total number of particle
m2.interp=ones(m2.nf,m2.np);    %pre-allocation for interp (1 for obtained by other methods and 0 for real)
m2.radius=zeros(m2.nf,m2.np);   %pre-allocation of radius
m2.r=zeros(m2.nf,m2.np,2);     %pre-allocation of position
m2.mass=zeros(m2.nf,m2.np);     %pre-allocation of mass
%assign all the values into m2
for i=1:m2.size
    m2.radius(data2(i,5)+1,data2(i,6)+1) = data2(i,4);
    m2.r(data2(i,5)+1,data2(i,6)+1,1)=data2(i,1);
    m2.r(data2(i,5)+1,data2(i,6)+1,2)=data2(i,2);
    m2.interp(data2(i,5)+1,data2(i,6)+1)=0;
    m2.mass(data2(i,5)+1,data2(i,6)+1) = data2(i,3);
end
data2=[];                                   %release the RAM of data2 occupied
'read data2 ok'

rmax=squeeze(max(max(m2.r,[],1),[],2));     %find the largest m.r2
rmin=squeeze(min(min(m2.r,[],1),[],2));     %find the smallest m.r2
L=(rmax-rmin);                              %set L
diameter=(L(1)*L(2)/m2.np).^0.5;             %diameter?????

% largest_diameter=max(m1.radius(m1.nf,:));            %set the smallest dr for particle

r1last=squeeze(m1.r(m1.nf,:,:));            %get the position of the last frame of data1
r2first=squeeze(m2.r(1,:,:));               %get the position of the first frame of data2

identifier=-1*ones(1,m1.np);                %set the identifier to be all -1
assigned_flag=zeros(1,m2.np);               %set the assignment flag to be all 0
datac=zeros(m1.size+m2.size,6);             %Pre-allocation of datac
l=1;                                        %initialize l to be the row number of the datac which keep adding up

for i=1:m1.np                               %for every particle in m1
    if m1.interp(m1.nf,i)==0                %this is a real particle
        drmin=diameter;                     %set the smallest dr
        for j=1:m2.np                       %for every particle in m2
            if assigned_flag(j)==0          %this particle has not been matched
                if m2.interp(1,j)==0            %if this is a real particle
                    dx=abs(r2first(j,1)-r1last(i,1));       %calculate the dx
                    if dx<drmin                 %if dx smaller than drmin
                        dy=abs(r2first(j,2)-r1last(i,2));   %calculate the dy
                        if dy<drmin                 %if dx and dy all less than drmin (Save computational power)
                            drmin=sqrt(dx*dx+dy*dy);%update drmin
                            identifier(i)=j;        %assign the j to the identifier(i) , points to the other particle in data2
                        end
                    end
                end
            end
        end
        if identifier(i)~=-1
            assigned_flag(identifier(i))=1;
        end
    end
    p=identifier(i);
    for f1=1:m1.nf                      %for every frame of this particle
        if m1.interp(f1,i)==0           %if this is a real particle
            datac(l,1)=m1.r(f1,i,1);    %set x into datac
            datac(l,2)=m1.r(f1,i,2);    %set y into datac
            datac(l,3)=m1.mass(f1,i);   %set mass into datac
            datac(l,4)=m1.radius(f1,i); %set radius into datac
            datac(l,5)=f1;              %set the frame number
            datac(l,6)=i;               %set the particle number
            l=l+1;
        end
    end
    if p~=-1                    %if it found one corresponding particle in r2last
        for f2=1:m2.nf                      %for every frame of this particle in data2
            if m2.interp(f2,p)==0           %if this is a real particle
                datac(l,1)=m2.r(f2,p,1);    %set x into datac
                datac(l,2)=m2.r(f2,p,2);    %set y pnto datac
                datac(l,3)=m2.mass(f2,p);   %set mass pnto datac
                datac(l,4)=m2.radius(f2,p); %set radpus pnto datac
                datac(l,5)=f2+m1.nf;        %set the frame number
                datac(l,6)=i;               %set the particle number
                l=l+1;
            end  
        end
    end
end

%add the un-paired particles of m2
for i=1:m2.np                               %for all the particles in m2         
    if assigned_flag(i)==0                  %if it is unpaired
        for f2=1:m2.nf                      %for every frame of this particle
            if m2.interp(f2,i)==0           %if this is a real iarticle
                datac(l,1)=m2.r(f2,i,1);    %set x into datac
                datac(l,2)=m2.r(f2,i,2);    %set y into datac
                datac(l,3)=m2.mass(f2,i);   %set mass into datac
                datac(l,4)=m2.radius(f2,i); %set radius into datac
                datac(l,5)=f2+m1.nf;        %set the frame number
                datac(l,6)=i+m1.np;         %set the particle number
                l=l+1;
            end
        end
    end
end

%slice it to real size
% datac=datac(1:find(datac(:,3)==0,1)-1,:);     %cut out the 0 terms
datac(:,5)=datac(:,5)-1;                    %-1 for all frame number (Since MATLAB starts from 1)
datac(:,6)=datac(:,6)-1;                    %-1 for all particle number (Since MATLAB starts from 1)

%Reset the sixth coloum to let it be a sequenced list
l=size(datac(:,6),1);                         %set the l be the #rows in datac
j=0;
last=datac(1,6);                        %initialize some parameters
for i=1:l                                   %for every row
    if last==datac(i,6)                     %if same as last row (the same particle)
        datac(i,6)=j;                       %make it to be j
    else                                    %this is another particle
        j=j+1;                              %add j
        last=datac(i,6);                    %store it to be last
        datac(i,6)=j;                       %make it to be j
    end
end
status=unix('rm 1.dat');
status=unix('rm 2.dat');
dlmwrite('c.dat', datac,'\t');
