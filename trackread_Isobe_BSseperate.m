%%This function is to read the Isobe data file using two files, one for
%%particle trajectories and the other for the first frame of the simulation
%%to give the radius and seperate the big and small particles
function [m,mbig,msmall] = trackread_Isobe_BSseperate(filename)

file_0=[filename num2str(0)];       %construct a file_0 name for the dat0 
data=dlmread(filename);             %read the posisiton file
data0=dlmread(file_0);              %read the radius file
1
dt=1;
[h, l] = size(data);                %get the size of data
m.d=2;
m.nframe=1000;                      %set the frame number as Parameter
m.n=4096;                           %set the particle number as Parameter
m.dt=dt;                            %initialize the m.fields
m.T = ones(m.nframe, 1) * 1;        %initialize the m.fields
m.Bulk=0;                           %initialize the m.fields

m.radius=zeros(m.nframe,m.n);       %initialize the m.fields
%set the radius of all particles
for i=1:m.nframe
    m.radius(i,:)=data0(:,3);
end

%sort
m.r=zeros(m.nframe,m.n,m.d);        %initialize the m.fields
m.interp=ones(m.nframe,m.n);        %initialize the m.fields

c=0;
for a=1:m.n;
    for b=1:m.nframe;
    m.r(b,a,1)=data(c+b,1);
    m.r(b,a,2)=data(c+b,2);
    m.interp(b,a)=0;
    end;
    c=c+m.nframe;
end;
1
%position
for p=1:m.n;
    last =0;
    ck=0;
    for f=1:m.nframe;
        if(m.interp(f,p)==1)
            if(ck==0)
                last=f;
                ck=1;
            end;
        elseif(ck==1)
            if(last==1)
                for tf=1:f-1;
                    m.r(tf,p,1)=m.r(f,p,1);
                    m.r(tf,p,2)=m.r(f,p,2);
                end;
            else
                dx = (m.r(f,p,1) - m.r(last - 1,p,1)) / (f - last + 1);
                dy = (m.r(f,p,2) - m.r(last - 1,p,2)) / (f - last + 1);
                for x = last:f-1;
                    m.r(x,p,1) = m.r(last - 1,p,1) + dx*(x - last + 1);
                    m.r(x,p,2) = m.r(last - 1,p,2) + dy*(x - last + 1);
                end;
            end;
            ck = 0;
            last = f;
        else
            last=f;
        end;
    end;
    if (ck == 1)
        for tf = last:m.nframe;
            m.r(tf,p,1) = m.r(last - 1,p,1);
            m.r(tf,p,2) = m.r(last - 1,p,2);
        end;
    end;
end;
m.r(:,:,3)=m.r(:,:,1)*0;
rmax=squeeze(max(max(m.r,[],1),[],2));
rmin=squeeze(min(min(m.r,[],1),[],2));
L=(rmax-rmin);
diameter=(L(1)*L(2)/m.n)^0.5;
m.r=m.r/diameter; % normalize to unit particle size
m.rmax=squeeze(max(max(m.r,[],1),[],2));
m.rmin=squeeze(min(min(m.r,[],1),[],2));
% flip origin from top-left to bottom left
m.r(:,:,2)=m.rmax(2)-m.r(:,:,2)+m.rmin(2);
% find stuck particles
% 0 for non-stuck, 1 for stuck
ng=10;
df=diameter*0.001;  % smaller than diameter*0.001
m.stuck=zeros(m.nframe,m.n);
for p=1:m.n;
    for f=m.nframe:-1:2;
        dx=m.r(f,p,1)-m.r(f-1,p,1);
        dy=m.r(f,p,2)-m.r(f-1,p,2);
        if(sqrt(dx*dx+dy*dy)<df)
            m.stuck(f,p)=1;
        else
            for tf=f+1:f+ng;
                if(tf>m.nframe)
                    break;
                else
                    m.stuck(tf,p)=0;
                end;
            end;
            break;
        end;
    end;
end;
%%Seperate the big and small
fields={'radius','r','interp','rmax','rmin','stuck'};
mbig=rmfield(m,fields);
msmall=rmfield(m,fields);

radius_small=min(data0(:,3));       %find the big and small radius of the particle
radius_big=max(data0(:,3));

mbig.radius=radius_big;
msmall.radius=radius_small;
big_identifier=(data0(:,3)==radius_big);

mbig.r=m.r(:,big_identifier,:);
mbig.radius=m.radius(:,big_identifier);
mbig.interp=m.interp(:,big_identifier);
mbig.n=sum(big_identifier);
mbig.stuck=m.stuck(:,big_identifier);

msmall.r=m.r(:,~big_identifier,:);
msmall.radius=m.radius(:,~big_identifier);
msmall.interp=m.interp(:,~big_identifier);
msmall.n=sum(~big_identifier);
msmall.stuck=m.stuck(:,~big_identifier);

mbig.rmax=squeeze(max(max(mbig.r,[],1),[],2)); 
mbig.rmin=squeeze(min(min(mbig.r,[],1),[],2)); 
msmall.rmax=squeeze(max(max(msmall.r,[],1),[],2)); 
msmall.rmin=squeeze(min(min(msmall.r,[],1),[],2)); 




