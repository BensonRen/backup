%function natom2 = showdisp(m,t,dt,varargin)
function [m,redprecentage] = showdisp_Ben(m,t,dt,varargin)
% 3D plot of show displacement vectors
% time is measured in frame
%
% function out = showdisp(m,t,dt,varargin)
% opt = getopt( struct('nlist',0, 'top',0, 'drmin', 0, 'Rnei', 0, 'n0',0, 'nseg',20, ...
%       'nowrap', 'noarg', 'chain', 'noarg'), varargin{:});

opt = getopt( struct('nlist',0, 'top', 0,'top_radius',0, 'drmin', 0, 'Rnei', 3, 'zmax',0, ...
       'n0',0, 'plotn0','noarg', 'nseg',50, ...
       'thinout','noarg', 'nowrap', 'noarg', 'chain', 'noarg', 'drmax',0, 'widthdr',0, 'MS', 10, 'LW',2, ... 
       'color',0, 'colory','noarg', 'colordfr',0, ...
       'save',0, 'drpath','noarg',  'accurate','noarg', 'drcut','noarg',  'line','noarg', ...
       'drawcnt',0, 'noclear','noarg', 'color0',0, ...
       'plotdim',0, ...
       'patch','noarg', 'patchn',4, 'sfinal',1, 'sinit',0, 'pub','noarg', 'movie','noarg', 'part',0, ...
       'arrow','noarg', ...              
       'dr2',1.2, 'link','noarg', 'drlink',0.4, ... 
       'dx',0, 'dy',0, 'dz',0, 'view', [], ...
        'xyt','noarg', 'jumponly','noarg', 'w',0.05, 'rsphere',0, 'nojoinsphere','noarg',  ...
       'sdrmax',0, 'cdrmax',0, 'cdrmin',0, 'head','noarg', 'showi','noarg', 'nlabel',[], ...
        'lscn','noarg','show_radius','noarg','real_cen','noarg'), varargin{:});

% Rnei: radius of neighborhood to show
% nseg: no. of time-points for each atom (>=2), equal no. of line seg plus 1 
% head: chain head or tail only

redprecentage=0; %initialize the redprecentage

% select particles %%%%%%%%%%%%%%%%%%%%%%%%%%%
w=opt.w;
if (opt.rsphere==0) opt.rsphere=w; end;
if (dt==0) dt = m.nframe-t; end;

[tmax natom d] = size(m.r);

t1=t;
t2=t+dt;
dt2=t2-t1+1;
%ds=m.s(t2)-m.s(t1)

T = mean(m.T(t1:t2));
r1=squeeze(m.r(t1,:,:));
r2=squeeze(m.r(t2,:,:));

if (opt.drpath) % consider dr along whole path
  Dt=t2-t1; Dtcnt=1;
  while (Dt>=1)
    if (Dt<10 || opt.accurate)
      Drseg = m.r(t1+Dt:1:t2,:,:)-m.r(t1:1:t2-Dt,:,:);
      Dt=Dt-1;
    else
      Drseg = m.r(t1+Dt:Dt:t2,:,:)-m.r(t1:Dt:t2-Dt,:,:); % skip some checks to speedup
      Dt=round(Dt/2.)-1;
    end;
    dr2seg=dot(Drseg,Drseg,3);
    dr2max(Dtcnt,:)=max(dr2seg,[],1);
    Dtcnt=Dtcnt+1;
  end;
  dr2=max(dr2max,[],1);
  dr=dr2.^0.5;
else % consider net displacements dr
  Dr=r2-r1;
  dr=dot(Dr,Dr,2).^0.5;
end;

if (opt.zmax>0)
  if (m.Bulk || m.SUBS>0)
    dr(r1(:,3)>opt.zmax)=0;
  else
    dr(abs(r1(:,3))>opt.zmax)=0;
  end;
end;

if (opt.drcut>0) % cut away atoms which moved at beginning and end 1/5 of the period
  Drseg = m.r(t1+1:t2,:,:)-m.r(t1:t2-1,:,:);
  drseg = dot(Drseg,Drseg,3).^0.5;
  dtcut=round((t2-t1)/5);
  drcut = max(drseg([1:dtcut end-dtcut:end],:),[],1);
  dr(drcut>opt.drmin)=0;
end;


if (opt.top>0) % show fastest "opt.top" percents
    drsort = sort(dr);
    opt.drmin = drsort( round(m.n*(1.-(opt.top/100.)) ));
end;

coloratom=ones(m.n,1);
nlist=dr>-1; % start with all atoms
disp(['size(nlist)', size(nlist)])
if (opt.nlist>0) % select atoms by number 
  coloratom=zeros(m.n,1);
  nlist = opt.nlist;
  coloratom(nlist)=1;
end;

  
if opt.drmin>0
  nlist = dr > opt.drmin; % all atoms with disp > opt.drmin
disp(size(nlist))
end;

if (opt.n0>0 && opt.Rnei>0) % neighborhood of atom opt.n0 at t1 
  w=w/4; % thinner cylinder
  r0=squeeze(m.r(t1,opt.n0,:));  disp(['r0 = ' num2str(r0')]);
  drn = squeeze(m.r(t1,:,:))-ones(natom,1)*r0';
  Rmax=ones(natom,1)*m.rmax;
  drn = mod(drn+Rmax/2,Rmax)-Rmax/2; % PBC in xyz
  drnnorm = dot(drn,drn,2).^0.5;
  neilist = drnnorm < opt.Rnei;
  if (size(neilist,1)~=size(nlist,1)) neilist=neilist'; end;
  nlist = nlist & neilist; % merge with AND
  opt.MS=opt.MS*4;
end;

if opt.drmax>0 % show only smaller than drmax
  nlist = dr < opt.drmax;
end;

if opt.head % show only head/tail of polymer chain
  ii = [1:m.n]';
  headlist = (mod(ii,m.len)<=1);
  size(nlist)
  size(headlist)
  nlist = and( nlist, headlist);
end;

% processing trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%

i = (1:m.n); atomi = i(nlist); 
%dr(nlist)' % DEBUG

if (opt.plotdim>0) % plot also r(opt.plotdim) vs t for at most 10 particles
  for (n=1:min(10,numel(atomi)))
    myplot('ctitle',str(atomi(n),0), 'lp', m.r(t1:t2,atomi(n),opt.plotdim));
    hold on;
  end;
  hold off;
  figure;
end;

if (opt.drawcnt>0) % for showdispmulti.m
  if (~isfield(m,'drawcnt') || numel(m.drawcnt)<m.n)  m.drawcnt = zeros(1,m.n); end;
  m.drawcnt(nlist) = m.drawcnt(nlist)+1;
  drawcntmax=max(m.drawcnt)
end;

% thinning-out/averaging data
rc0 = m.r(t1:t2,nlist,:);
if (opt.n0) n0=sum(nlist(1:opt.n0)); end;
dr = dr(nlist);
if (opt.nseg*2>dt2) 
  rc=rc0;
elseif opt.thinout % thin out path positions
    seglist=0:opt.nseg;
    rc(:,:,:)=rc0(int32(seglist*dt/opt.nseg)+1,:,:); 
else % sub-average path positions
  ds=dt2/opt.nseg;
  for tcnt=1:opt.nseg
    rc(tcnt,:,:)=mean(rc0(uint32((tcnt-1)*ds)+1:uint32(tcnt*ds),:,:),1);
  end;
end;
if (size(rc,3)==1) rc=reshape(rc,size(rc,1),1,3); end; % 1 dimension missed due to 1 ptcle selected

if (isfield(m,'interp'))
  tm=round(mean(t1,t2));
  interp1 = min(m.interp(t1:tm,nlist),[],1);
  interp2 = min(m.interp(tm:t2,nlist),[],1);
  interp = (interp1==1) | (interp2==1); % either first or second half data missing
end;

if (isfield(m,'stuck'))
  tm=round(mean(t1,t2));
  stuck1 = min(m.stuck(t1:t2,nlist),[],1);
  stuck = (stuck1==1) ; 
end;

if (opt.n0 & opt.link) % show only ptcles linked (i.e. closeby) to n0 
  linked = findlinked(rc,n0,opt.drlink,[]); 
  nlinked=sum(linked) %DEBUG
  m.nlinked=nlinked;
  atomi;
  linked;
  atomlist=atomi(linked==1) %DEBUG: show selected atoms
end;

[tmax2 natom2 d] = size(rc);
disp(['[ tmax2 natom2 d ] = ' num2str(size(rc)) ]);
m.natom2 = natom2;

if (opt.dx~=0) rc(:,:,1)=rc(:,:,1)+opt.dx; end; % shifting
if (opt.dy~=0) rc(:,:,2)=rc(:,:,2)+opt.dy; end;
if (opt.dz~=0) rc(:,:,3)=rc(:,:,3)+opt.dz; end;

if (natom2>0 && opt.n0==0) % pbc wrapping 
  if (~opt.nowrap) 
    if (m.Bulk==1)
      rc = wrappbc(rc,m.rmax,'dim',[1 2 3]); 
    else
      rc = wrappbc(rc,m.rmax,'dim',[1 2]); 
    end;
  end; 
elseif (opt.n0) 
  rc = wrapr0(rc,m.rmax,r0); % wrap in x,y directions
end; 

if (opt.n0 && opt.plotn0) % plot also x,y,z vs t
  for (n=n0:n0)
    for (dir=1:3)
      myplot('ctitle',str(atomi(n),0), 'lp', rc(:,n,dir));
      hold on;
    end;
  end;
  hold off;
  figure;
end;

rcmax = squeeze(max(max(rc,[],1),[],2));
rcmin = squeeze(min(min(rc,[],1),[],2));
L=max(rcmax-rcmin);
if (L>40) opt.MS=opt.MS/2; opt.LW=opt.LW/2; end; % small markers and thin lines

% draw path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (opt.noclear) 
  hold on; 
else
  SS = get(0,'screensize');
  if (opt.lscn)
    %set(0,'defaultfigureposition',[SS(3)-1100-50 SS(4)-600-50 1500 1100]);
    set(0,'defaultfigureposition',[0 1920-1080 1920 960]); % left full screen of apricot [x y w h]
  end;
  
  if (opt.patch) % full 3D lighting (slow)
    figure('Color','white'); 
    %axis off; 
    set(gca,'Xtick',0:-1); 
    set(gca,'Ytick',0:-1); 
    set(gca,'Ztick',0:-1); 
    box on;
    set(gca,'DataAspectRatio',[1 1 1],'Projection','perspective');
    set(gca,'DataAspectRatio',[1 1 1],'Projection','orthographic');
    campos([100 0 0]);
  else  % quick 3D
    clf;  
    axes; 
    set(gca,'DataAspectRatio',[1 1 1],'Projection','orthographic');
    set(gca,'Position',[0.02 0.02 0.98 0.98]);
  end;
end;
if (natom2==0)  return; end; % no particle selected

rotate3d on; % reset(zoom);
hold on;
set(gca,'Clipping','off');
%cameratoolbar; 
ncolor=1024;
ccut1=100; ccut2=100;  % cut away some dark color in jet palette
if (opt.color0>0) ccut1=round(opt.color0*(ncolor-1)); ccut2=0; end;
cmatrix=jet(ncolor+ccut1+ccut2);
cmatrix=cmatrix(ccut1:ccut1+ncolor-1,:);
%cadd2=350; cmatrix(end-cadd2:end-1,:)=ones(cadd2,1)*cmatrix(end,:); % add more red to for strint.jpg
graymap=colormap(gray(ncolor));
cmap=colormap(flipud(cmatrix));

disp('create graphic object G');
idr2=[];
Gc = [];
Gs = [];
if 1
  % show trajectories
  if (opt.n0 && n0>0) scatter3(rc(1,n0,1),rc(1,n0,2),rc(1,n0,3),1.5*opt.MS,'r','filled'); end;
  for (n=1:natom2)
    if (~opt.n0 || ~opt.link || linked(n))
      kmid=0;
      if (opt.sinit>0)  kmid=1+round((tmax2-1)*opt.sinit); end;
      if (opt.sfinal<1) kmid=1+round((tmax2-1)*opt.sfinal); end;
      cmid = 1+round(kmid/(tmax2)*(ncolor-1));
      for k = 1:tmax2-1;
        Rn = squeeze(rc(k:k+1,n,:));
        DRn=Rn(2,:)-Rn(1,:);
        dRn=norm(Rn(2,:)-Rn(1,:)); dRnk(k)=dRn;
        if (opt.colory || opt.xyt) 
          c = 1+round( mod(Rn(1,2)/m.rmax(2),1)*(ncolor-1) ); % color depending on y co-ord
        elseif (opt.colordfr>0) 
          c = 1+round( mod(t+k,opt.colordfr)/opt.colordfr*(ncolor-1) ); % color depending on t, cycled with period colordfr
        elseif (opt.drawcnt>0) % for showdispmulti.m
          c = 1+round(min(ncolor,max(1,(ncolor-1)*(1-m.drawcnt(atomi(n))/opt.drawcnt))));
        elseif (opt.color>0) 
          c=1+round(opt.color*(ncolor-1)); 
        else
          c = 1+round((k+0.5)/(tmax2)*(ncolor-1)); % 0.5 for cylinder between k & k+1
        end;
        if (opt.xyt) Rn(1,2)=k*0.001; Rn(2,2)=(k+1)*0.001; end; % plot trajectories in x-y-t coords (rather than x-y-z)
        if (~opt.nowrap) % if wrapped across pbc, don't show 
          if(abs(Rn(1,1)-Rn(2,1))>m.rmax(1)/2||abs(Rn(1,2)-Rn(2,2))>m.rmax(2)/2||abs(Rn(1,3)-Rn(2,3))>m.rmax(3)/2) c=0; end;
        end;
        if (c>0 && (~opt.xyt || dRn>opt.drmin) && (~opt.jumponly || dRn>opt.drmin))
          if (coloratom(atomi(n)))
            cm = cmap(c,:);  % color
          else
            cm = graymap(c,:); % gray
          end;
          %if (opt.xyt&&dRn<opt.drmin) cm = graymap(round(0.9*ncolor),:); end;
          if (opt.arrow)
            rx=dot(DRn,[1 0 0])/dRn; rz=dot(DRn,[0 0 1])/dRn; ry=(rx*rx+rz*rz)^0.5;
            arrow(Rn(1,:), Rn(2,:),'EdgeColor',cm, 'FaceColor',cm, 'Length',10*ry, 'TipAngle',min(60,20/ry) , 'CrossDir',cross([0 1 0],DRn), 'Width',2 );
          elseif (~opt.patch)
            color(cm);
            line(Rn(:,1), Rn(:,2), Rn(:,3), 'color', cm, 'LineWidth',opt.LW);
          else
            %G = gobjadd(G, cylinderpatch(w, 3, 2, Rn(1,:), Rn(2,:), c)); 
            if opt.movie && opt.part>0
              RA=Rn(1,:); RB=Rn(2,:); RBA=RB-RA;
              npart=7; part=opt.part;
              cA=max(0,(part-4)/(npart-3));
              cB=min(1,(part)/(npart-3));
              Ra=RA+cA*RBA;
              Rb=RA+cB*RBA;
              Gc = gobjadd(Gc, cylinderpatch(w, opt.patchn, 2, Ra , Rb, c)); 
              Gs = gobjadd(Gs, spherepatch(w, opt.patchn, Ra', c));
              Gs = gobjadd(Gs, spherepatch(w, opt.patchn, Rb', c));
            else
              Gc = gobjadd(Gc, cylinderpatch(w, opt.patchn, 2, Rn(1,:), Rn(2,:), c));
            end
            if (~opt.nojoinsphere && ~opt.movie) Gs = gobjadd(Gs, spherepatch(w, opt.patchn, Rn(1,:)', c)); end;
          end;
        end;
      end;
      if (c>0 && opt.rsphere>0 && ~opt.xyt) % draw/add graphic object
        if (~opt.patch)
          if (isfield(m,'interp') && interp(n))
            scatter3(rc(1,n,1),rc(1,n,2),rc(1,n,3),opt.MS*4,'k','filled','o'); % initial position
          elseif (isfield(m,'stuck') && stuck(n))
            scatter3(rc(1,n,1),rc(1,n,2),rc(1,n,3),opt.MS*4,'b','filled','o'); % initial position
          else
            %color(cm);
            scatter3(rc(1,n,1),rc(1,n,2),rc(1,n,3),opt.MS,'r','filled','o'); % initial position
            %scatter3(rc(end,n,1),rc(end,n,2),rc(end,n,3),opt.MS,'b','filled','<'); % final position
          end;
        else
          %Gs=gobjadd(Gs, spherepatch(opt.rsphere   , 1*opt.patchn, rc(1,n,:), 1)); % initial position
          Gs=gobjadd(Gs, spherepatch(opt.rsphere, 1*opt.patchn, rc(1,n,:), c)); % initial position
          if (kmid>0) Gs=gobjadd(Gs, spherepatch(opt.rsphere*0.99   , 1*opt.patchn, rc(kmid,n,:), cmid)); end;% intermediate
        end;
      end;
      if (opt.dr2>0 && dr(n)>opt.dr2 && ~opt.xyt)
        scatter3(rc(1,n,1),rc(1,n,2),rc(1,n,3),2*opt.MS,'k','filled'); % big hop
        idr2(end+1)=atomi(n);
      end;
      if (opt.showi) 
        text(double(Rn(end,1))-0.3,double(Rn(end,2)),double(Rn(end,3)),str(atomi(n),0));
      end;
      if (numel(opt.nlabel)>=n) 
        text(double(Rn(end,1))-0.3,double(Rn(end,2)),double(Rn(end,3)),str(opt.nlabel(n),0));
      end;
    end;
  end;
end;
m.dRnmax=max(dRnk);

if (opt.dr2>0)  out=idr2'; save 'dispi2.txt' out -ascii; end;
%G = gobjadd(G, cylinderpatch(0.2, 6, 2, [0,0,0], [2,0,0], cmap(1,:))); % draw color bar
%G = gobjadd(G, cylinderpatch(0.2, 6, 2, [2,0,0], [4,0,0], cmap(ncolor,:))); 
  
disp('calling patch'); 
Gprop={'FaceColor','interp', 'CDataMapping','direct', 'FaceLighting','phong', 'FaceAlpha','interp','AlphaDataMapping','none', ... %'FaceAlpha',1., ...
       'AmbientStrength',0.5, 'DiffuseStrength',0.5, 'SpecularStrength', 0.5, 'SpecularExponent',15, 'EdgeColor','flat'};
if (numel(Gc)>0 && opt.w>0) 
  Gc.A = 1+0.*(Gc.C<=1);  % solid
  patch('Vertices', Gc.V, 'Faces', Gc.F, 'FaceVertexCData',Gc.C, 'FaceVertexAlphaData',Gc.A, Gprop{:});
end;
if (numel(Gs)>0 && opt.rsphere>0) 
  Gs.A=1+0.*(Gs.C<=1); % solid
  patch('Vertices', Gs.V, 'Faces', Gs.F, 'FaceVertexCData',Gs.C, 'FaceVertexAlphaData',Gs.A, Gprop{:}); 
end;
%axis tight; 
%set(gca,'visible','off'); set(gcf,'color',[1 1 1]);

hold off;
%view([0 0]); 
view(2);
d=m.rmax(1)*0.05; % margin width
if (opt.n0==0 && opt.xyt==0 && m.Bulk==1) 
  xlim([0-d m.rmax(1)+d]); ylim([0-d m.rmax(2)+d]); zlim([0-d m.rmax(3)+d]); 
  %xlim([0-d m.rmax(1)+d]); ylim([0-d m.rmax(1)+d]); zlim([-m.rmax(1)/2-d m.rmax(1)/2+d]); 
else
  %l=[-opt.Rnei,opt.Rnei]; xlim(m.r(t1,opt.n0,1)+l);  ylim(m.r(t1,opt.n0,2)+l);  zlim(m.r(t1,opt.n0,3)+l); 
end;
if (opt.zmax>0) && ~opt.noclear
  xlim([0-d m.rmax(1)+d]); ylim([0-d m.rmax(1)+d]); zlim([-opt.zmax-d opt.zmax+d]); 
end;
  
%if (opt.xyt) view([0.1 0.1 0.1]); end;
%xlim([0 m.rmax(1)]); ylim([0 m.rmax(1)]); zlim([-15 15]); 
%grid on;
if (numel(opt.view)>0) view(opt.view); end;
if (opt.patch) 
  if opt.pub
    shading interp; 
    lighting gouraud; % for curved surfaces
    camlight('headlight');
    camlight('right');
    camlight(-30, -80); camlight(20, 40); camlight(25, -40); camlight(-110, 0);
    material([0.43 0.35 0.1]); 
  else
    shading interp; 
    lighting gouraud; % for curved surfaces
    %light('Position',[-1 -1 0], 'style','local'); 
    %camlight('headlight');
    if (~opt.movie) camlight(20, 40); camlight(25, -40); camlight(-110, 0); end;
    %material([0.26 0.26 0.25]); % str.jpg
    material([0.3 0.3 1]); % strrep.jpg
  end;
end;
%if (SS(4)>1080) zoom(1.2); else zoom(1.2); end;
%axis vis3d; 
%colorbar;
grid off;
rotate3d on; 
%if (opt.patch) axis off; end;
drawnow;






if (opt.show_radius>0)
    %hold on
    
    radius = transpose(me.radius(t,:));
    
    par = readtable('../density/parameter_radius.txt','Delimiter','\t');
    %radius_ratio = par.Value(1);
    %dradius = (m.radius_big-m.radius_small*radius_ratio)/(radius_ratio-1)

    %zero = (radius == 0);
    %radius(~zero) = (radius(~zero) + dradius) .* m.radius_small ./ (m.radius_small + dradius);
    
    %radius_normalized = radius ./ m.radius_small;
    radius_small_real = par.Value(2);
    radius_normalized = radius ./ radius_small_real;
   
    
    if (opt.real_cen > 0)
        centres = squeeze(m.r(t,:,1:2));
    else
        %rc0 = m.r(t1:t2,:,:);
        ds = dt2 / opt.nseg;
        if uint32(ds) < 1
            ds = 1;
        end;
        %rc(1,:,:) = mean(rc0(1: uint32(ds), :, :), 1);
        for tcnt=1:opt.nseg
            rc(tcnt,:,:)=mean(rc0(uint32((tcnt-1)*ds)+1:uint32(tcnt*ds),:,:),1);
        end;
        %centres = squeeze(rc(1,:,1:2));
        centres_all = permute(rc,[2,3,1]);
        %centres = centres(:,:,1);
    end
    
    if (opt.top_radius>0)
        %rc0 = m.r(t1:t2,:,:);
        ds=dt2/opt.nseg;
        %for tcnt=1:opt.nseg
        %    rc(tcnt,:,:)=mean(rc0(uint32((tcnt-1)*ds)+1:uint32(tcnt*ds),:,:),1);
        %end;
        Drc = diff(rc,1,1);
        drc2 = dot(Drc,Drc,3);
        %KE_av = sum(drc2,1)/(opt.nseg-1);
        KE_av = squeeze(sum(drc2,1));
        %Drc0 = diff(rc0,1,1);
        %drc0_2 = dot(Drc0,Drc0,3);
        %KE_av = squeeze(sum(drc0_2,1));
        dr = KE_av .^ 0.5;
        drsort = sort(dr);
        seg_top = 10;

        magnification = 1 / 2;

        drmin_radius = drsort( round(m.n*(1.-(opt.top_radius/100.)) ));
        nlist_radius = (transpose(dr) <= drmin_radius);
        centres = centres_all(:,1:2,1);
        viscircles(centres(nlist_radius,:),magnification*radius_normalized(nlist_radius),'linewidth',0.5,'linestyle','-');
        hold on

        %cmap = jet(seg_top*100/opt.top_radius);    %full range
        %cmap = parula(seg_top);
        cmap = winter(seg_top);
        %cmap = [0,0,1;0,1,1;0,1,0];
        for i = 1:seg_top
            drmin_radius = drsort( round(m.n*(1.-(i*(opt.top_radius/seg_top)/100.)) ));
            drmax_radius = drsort( round(m.n*(1.-((i-1)*(opt.top_radius/seg_top)/100.)) ));
            nlist_radius = (transpose(dr) > drmin_radius) & (transpose(dr) <= drmax_radius);
            %centres_radius = centres(nlist_radius,:);
            %circle(centres_radius(:,1),centres_radius(:,2),magnification*radius_normalized(nlist_radius),cmap(i,:));
            viscircles(centres(nlist_radius,:),magnification*radius_normalized(nlist_radius),'EdgeColor',cmap(i,:));
        end;
        
        axis('equal')
    else
        nseg_rc = round(dt/10)
        if nseg_rc < 2
            nseg_rc = 2
        end;
        %rc0 = m.r(t1:t2,:,:);
        ds=dt2/nseg_rc;
        for tcnt=1:nseg_rc
            rc_dr(tcnt,:,:)=mean(rc0(uint32((tcnt-1)*ds)+1:uint32(tcnt*ds),:,:),1);
        end;
        t1 = 1;
        t2 = nseg_rc;
        dt2=t2-t1+1;
        Dt = t2 - t1;
        Dtcnt = 1;
        while (Dt>=1)
            Drseg = rc_dr(t1+Dt:1:t2,:,:)-rc_dr(t1:1:t2-Dt,:,:);
            Dt=Dt-1;
            dr2seg=dot(Drseg,Drseg,3);
            dr2max(Dtcnt,:)=max(dr2seg,[],1);
            Dtcnt=Dtcnt+1;
        end;
        dr2=max(dr2max,[],1);
        dr=dr2.^0.5;

        drmin_radius = 0.6;

        nlist_radius = transpose(dr) >= drmin_radius;
        
        magnification = 1 / 2;
        m.nseg = opt.nseg;
        
        for frame = (1:opt.nseg) 
            centres = centres_all(:,1:2,frame);
            centres_string = centres(nlist_radius,:);
            centres_nonstring = centres(~nlist_radius,:);
            radius_string = radius_normalized(nlist_radius,:);
            radius_nonstring = radius_normalized(~nlist_radius,:);
            redprecentage=length(centres_string)/(length(centres_all));
            hold on     
            %c1 = circle(centres_string(:,1),centres_string(:,2),magnification*radius_string,'b');
            %c1 = viscircles(centres(nlist_radius,:),magnification*radius_normalized(nlist_radius),'linewidth',0.5,'linestyle','-','EdgeColor','r');
            
            %Add a percentage of the red particles into this
            text(5,2,['Red particle percentage:',num2str(redprecentage)]);
            
            for i = (1:length(centres_string))
                r_temp = radius_string(i)*magnification;
                c_temp = centres_string(i,:);
                %pos_temp = [c_temp-r_temp 2*r_temp 2*r_temp];
                %rec = rectangle('Position', pos_temp, 'Curvature', [1 1], 'FaceColor', 'r', 'Edgecolor', 'r');
                cir = filledCircle(c_temp, r_temp, 500, 'w'); %This is to fill the circle
                %color the circle 
                alpha(cir, 0.001);
                set(cir, 'EdgeColor', 'r');
                hold on
            end

            hold on
            %c2 = viscircles(centres(~nlist_radius,:),magnification*radius_normalized(~nlist_radius),'linewidth',0.5,'linestyle','-', 'EdgeColor', [0.5 0.5 0.5]);
            for i = (1:length(centres_nonstring))
                r_temp = radius_nonstring(i)*magnification;
                c_temp = centres_nonstring(i,:);
                %pos_temp = [c_temp-r_temp 2*r_temp 2*r_temp];
                %rec = rectangle('Position', pos_temp, 'Curvature', [1 1], 'FaceColor', [0.5 0.5 0.5], 'Edgecolor', [0.5 0.5 0.5]);
                cir = filledCircle(c_temp, r_temp, 500, 'b');
                alpha(cir, 0.05);
                set(cir, 'EdgeColor', 'b');
                hold on
            end
            axis('equal')
            if (m.gif == 0 && frame == 1)
                break
            end
            hold on;
            progress = text(min(centres(:,1,1)),min(centres(:,2,1)),[num2str(frame),'/',num2str(opt.nseg)]);
            %grid off;
            drawnow;
            %pause(0.1)

            fff = getframe(1);
            m.im{frame} = frame2im(fff);

            delete(progress);
            delete(c2)
            delete(c1)
            %for k = 1:length(c1)
            %    delete(c1{k})
            %end
        end
    end;

    
    %hold on
    %plot(centres(:,1),centres(:,2),'r.');
    hold off

end;

