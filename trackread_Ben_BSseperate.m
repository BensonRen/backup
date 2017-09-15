function [m,mbig,msmall] = trackread_Ben_BSseperate(filename)
% Create m from trackoutput.dat
% m.interp=0 for real position, m.interp=1 for the ones
% obtained by some techniques
dt=1;
data = dlmread(filename);

%filter extremely short path    %Ping
%%%%par_label = data(:,6);
%%%%appear_label = [];
%%%%label = 0;
%%%%data(end,end)
%%%%frame_max = max(data(:,5))
%%%%for i = (0:data(end,end))
%%%%    appear = sum(par_label==i);
%%%%    par_label = par_label(appear+1:end);
%%%%    %min(par_label)
%%%%    appear_label = [appear_label; ones(appear,1) * appear];
%%%%    if appear >= frame_max*0.1
%%%%        temp = [temp; ones(appear,1) * label];
%%%%        label = label + 1;
%%%%    end;
%%%%end;
%%%%data = data(appear_label>=frame_max*0.1,:);
%%%%data(:,6) = temp;
%%%%data(end,end)

[h, l] = size(data);

m.d=2;                      %dimension
%%%%m.nframe=data(h,l-1)+1;     %total number of frame of the last particle
m.nframe=max(data(:,l-1))+1;
m.n=data(h,l)+1;            %total number of particle
m.dt=dt;
m.T = ones(m.nframe, 1) * 1;
m.Bulk=0;
%sort
%Pre-allocation for the three big matrix
m.radius = zeros(m.nframe,m.n);
m.r=zeros(m.nframe,m.n,m.d);    %position of particles in different frames
m.interp=ones(m.nframe,m.n);

%put particle position, radius into m.radius and m.r
%set interp=0 for all the particles which have real positions by tracking (others remains one)
for i=1:h;
    %%%%m.radius(data(i,l-1)+1,data(i,l)+1) = sqrt(data(i,4));
    m.radius(data(i,l-1)+1,data(i,l)+1) = data(i,4);
    m.r(data(i,l-1)+1,data(i,l)+1,1) = data(i,1);
    m.r(data(i,l-1)+1,data(i,l)+1,2) = data(i,2);
    m.interp(data(i,l-1)+1,data(i,l)+1) = 0;
end;
'ok1'

%averaging the radii to avoid out-focus     %Ping
%for every particle
for i = 1:m.n
    %radius_temp is the radius in all frame of that particle
    radius_temp = m.radius(:,i);
    %create logical sequence for those 0s (Missing of particle)
    zero = ( radius_temp == 0 );
    %calculate the mean of the non-zero radius
    mean_non_zero = mean(radius_temp(~zero));
    %assign all the non-zero radius with same mean value
    radius_temp(~zero) = mean_non_zero;
    %assign back to teh m struct
    m.radius(:,i) = radius_temp;
    %%%%radius_all(i) = mean_non_zero;
end;

%position
for p=1:m.n;    %particle label, go through all the particles 
    last =0;    
    ck=0;
    for f=1:m.nframe;   %frame label, go through all the frame3
        if(m.interp(f,p)==1)    %if this particle is not obtained by IDL
            if(ck==0)
                last=f;     %update the 'last' value
                ck=1;
            end;
        elseif(ck==1)
            if(last==1)     
                %if there is no known value ahead(last=1), %use the next value to replace
                for tf=1:f-1;
                    m.r(tf,p,1)=m.r(f,p,1);
                    m.r(tf,p,2)=m.r(f,p,2);
                end;
            else
                %compute on average, dx, dy for the particle for 1 frame of
                %time
                dx = (m.r(f,p,1) - m.r(last - 1,p,1)) / (f - last + 1);
                dy = (m.r(f,p,2) - m.r(last - 1,p,2)) / (f - last + 1);
                %assume the particle moves the same length at each time
                %interval, calculate all the middle values of the particle
                for x = last:f-1;
                    m.r(x,p,1) = m.r(last - 1,p,1) + dx*(x - last + 1);
                    m.r(x,p,2) = m.r(last - 1,p,2) + dy*(x - last + 1);
                end;
            end;
            ck = 0;
            last = f;       %update the 'last' value
        else
            last=f;     %update the 'last' value
        end;
    end;
    if (ck == 1)
        %set all the 0s in m.r to be the same as the last non-zero value of m.r
        for tf = last:m.nframe;
            m.r(tf,p,1) = m.r(last - 1,p,1);
            m.r(tf,p,2) = m.r(last - 1,p,2);
        end;
    end;
end;
'ok2'

%make one more dimension for m.r and make it all 0s
m.r(:,:,3)=m.r(:,:,1)*0;
%%%%rmax=squeeze(max(max(m.r,[],1),[],2));
%%%%rmin=squeeze(min(min(m.r,[],1),[],2));
%%%%L=(rmax-rmin);
%%%%diameter=(L(1)*L(2)/m.n)^0.5;

%%=This is for plotting the histogram and fit it with a guassian function
%%as well as the mid point dividing the big and small particles
%===========================================================================
radius_first = m.radius(1,:);
radius_origin=radius_first;
%h = histogram(radius_first(within),1000);
figure;
h = histogram(radius_first(radius_first~=0),100);
%h = histogram(radius_all(within),1000);
bins = transpose((h.BinEdges(1:end-1)+h.BinEdges(2:end))./2);
f = fit(bins,transpose(h.Values),'gauss2');
hold on
plot(f);
if f.b1 > f.b2
    m.radius_big = f.b1
    m.radius_small = f.b2
else
    m.radius_small = f.b1
    m.radius_big = f.b2
end;

xg = (m.radius_small:0.01:m.radius_big);
yg = f.a1.*exp(-((xg-f.b1)./f.c1).^2) + f.a2.*exp(-((xg-f.b2)./f.c2).^2);
m.boundary = xg(find(yg == min(yg)))
hold on
plot([m.boundary, m.boundary], [min(yg), max(yg)], 'linewidth', 3);
drawnow;
screen2png('radius_distribution')
%==========================================================================

    
%%%%radius_sort = sort(radius_first);      %for isobe data
%%%%m.radius_small = radius_sort(1);       %for isobe data
%%%%m.radius_big = radius_sort(end);       %for isobe data
%%%%magnification = sqrt(2);                %for isobe data

%%%%par = readtable('parameter_radius.txt','Delimiter','\t');
par = readtable('../density/parameter_radius.txt','Delimiter','\t');       %read the parameter_radius file 
m.radius_ratio = par.Value(7)       %read the radius ratio
radius_small_real = par.Value(2);   %read the small particle radius
radius_big_real = par.Value(2) * m.radius_ratio;    %compute the big particle raidus as small*ratio
radius_add_small = radius_small_real-m.radius_small     
%===========???????============
%compute the difference in small particle radius between the first frame and the previously computed value? Why?
radius_add_big = radius_big_real-m.radius_big
%add the big_add and small_add to the big and small particles seperately
m.radius(m.radius~=0 & m.radius>m.boundary) = m.radius(m.radius~=0 & m.radius>m.boundary) + radius_add_big;
m.radius(m.radius~=0 & m.radius<=m.boundary) = m.radius(m.radius~=0 & m.radius<=m.boundary) + radius_add_small;

%Open a new figue and plot again the above histogram
%However, using the newly computed radius
%which has been added small adjustment on it
%================================================
figure;     
radius_first = m.radius(1,:);
h = histogram(radius_first(radius_first~=0),100);
screen2png('radius_distribution_after')
%================================================


%Compute the density in the region 'within'
%================================================
x = m.r(1,:,1);     %assign x postitions of the particles in the first frame
y = m.r(1,:,2);     %assign y postitions of the particles in the first frame
%define the region 'within'
within = (x>=par.Value(3)) & (x<=par.Value(4)) & (y>=par.Value(5)) & (y<=par.Value(6));
%calculate the density inside the 'within' area
density = sum(radius_first(within) .^ 2 *pi) / (par.Value(4)-par.Value(3)) / (par.Value(6)-par.Value(5))
%================================================

%%%%diameter = 2 * m.radius_small;
%set the diameter
diameter = 2 * radius_small_real;
%diameter = 2 * radius_small * magnification;
m.r=m.r/diameter; % normalize to unit particle size
%remove the singleton dimension
m.rmax=squeeze(max(max(m.r,[],1),[],2));    %find the largest x and y
m.rmin=squeeze(min(min(m.r,[],1),[],2));    %find the smallerst x and y
% flip origin from top-left to bottom left
m.r(:,:,2)=m.rmax(2)-m.r(:,:,2)+m.rmin(2);
% find stuck particles
% 0 for non-stuck, 1 for stuck
ng=10;
%df=diameter*0.001;  % smaller than diameter*0.001
df=0.001;  % smaller than 0.001
m.stuck=zeros(m.nframe,m.n);
for p=1:m.n;        %for every particle
    for f=m.nframe:-1:2;    %from the last frame to the second frame
        dx=m.r(f,p,1)-m.r(f-1,p,1);     %calculate dx (move in x)
        dy=m.r(f,p,2)-m.r(f-1,p,2);     %calculate dy (move in y)
        if(sqrt(dx*dx+dy*dy)<df)    %if it moves less than df
            m.stuck(f,p)=1;     %Then it is stuck
        else    %if it did not stuck
            for tf=f+1:f+ng;    %for the next ng frames,make them unstuck
                if(tf>m.nframe)
                    break;
                else
                    m.stuck(tf,p)=0;
                end;
            end;
            break;
        end;
    end;

end

%make a mbig and msmall for the same purpose, dividing the m into two
%catagories
fields={'radius','r','interp','rmax','rmin','stuck'};
mbig=rmfield(m,fields);
msmall=rmfield(m,fields);

big_identifier=(radius_origin>=m.boundary);

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
end



