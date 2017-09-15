%function showdispall(m,mm,dt,savename, varargin)
function showdispall_Ben(m,dt,savename, varargin)
[opt vararg] = getopt( struct('gif','noarg'), varargin{:});
clf;
nframe=size(m.r,1);
dir = [ 'disp' num2str(savename) '-' num2str(dt) ];
mkdir(dir);

dfr=dt;
if (isfield(m,'Ckmax') && m.Ckmax>1)
  dfr=m.Ckmax;
end;

SS = get(0,'screensize');
%set(0,'defaultfigureposition',[SS(3)-1000-50 SS(4)-900-50 800 600]);
%set(0,'defaultfigureposition',[SS(3)-1000-50 SS(4)-900-50 1400 900]);
set(0,'defaultfigureposition',[100 50 1800 1200]);

string_redpercentage=zeros(floor((nframe)/dfr)+1);
redpre_index=1;

for t=1:dfr:min(10*dfr,nframe-1)
%for t=1:dfr:nframe-1
  m.gif = opt.gif;
  if (opt.gif > 0)
      fname = [dir '/disp' num2str(savename) '-' num2str(dt) '-' num2str(t,'%05.0f') '.gif']
  else
      fname = [dir '/disp' num2str(savename) '-' num2str(dt) '-' num2str(t,'%05.0f') '.png']
  end
  if (~exist(fname,'file'))
    close all; 
    


    dt1=min(dt,nframe-t-1);
    [m2,redprecentage] = showdisp_Ben(m,t,dt1, vararg{:});
    %natom2 = showdisp_Ben(m,t,dt, 'nseg',20000, 'drmin',0.3, 'line', 'nowrap', 'colory', 'drseg');
    %natom2 = showdisp_Ben(m,t,dt, 'nseg',20000, 'drmin',0.3, 'line', 'nowrap', 'dr2',1.2 ,'colory', 'drseg');
    %natom2 = showdisp_Ben(m,t,dt, 'nseg',20000, 'drmin',0.3, 'line', 'dr2',1.2 ,'jumponly');
    %natom2 = showdisp_Ben(m,t,dt, 'nseg',200, 'drmin',0.6, 'line', 'dr2',10 ,'drmin',0.6, vararg{:});
    %natom2 = showdisp_Ben(m,t,dt, 'nseg',200, 'drmin',0.6, 'line', 'dr2',10 ,'drmin',0.6, 'dx',rand*m.rmax(1),'dy',rand*m.rmax(1), vararg{:});
    %m2 = showdisp_Ben(m,t,dt, 'nseg',200, 'drmin',0.6, 'line', 'dr2',10 ,'drmin',0.6, vararg{:});
    %m2 = showdisp_Ben(m,t,dt, 'nseg',200, 'drmin',0.6, 'line', 'dr2',10 ,'drmin',0.6, 'patch', 'patchn',10,'pub', 'w',0.1, 'colory','drpath',vararg{:});
    %m2 = showdisp_Ben(m,t,dt, 'nseg',2000, 'drmin',0.6, 'line', 'dr2',10 ,'drmin',0.6, 'patch', 'patchn',10,'pub', 'w',0.1, 'drpath', 'colordfr',dt, vararg{:});
    
    if (opt.gif > 0)
        gif_name = fname
        for idx = 1:m2.nseg                                                   
            [A, map] = rgb2ind(m2.im{idx}, 256);                                  
            if idx == 1                                                        
                imwrite(A, map, gif_name, 'gif', 'LoopCount', Inf, 'DelayTime', 1)    ;
            else                                                               
                imwrite(A, map, gif_name, 'WriteMode', 'append', 'DelayTime', 1);
            end                                                                
        end 
    else
        grid off;
        drawnow;
        pause(0.1);

        %theta=30+t*0.005*10; view([theta 5]); % rotate
        %view([0 1 0]);
        %view([0 0 1]);

        screen2png(fname);
        
        string_redpercentage(redpre_index)=redprecentage;
        redpre_index=redpre_index+1;
    end
  end;
      
end;

if (sum(cellfun(@(x)isequal(x,['show_radius']),varargin))>0)
    cd(dir);
    dlmwrite('percentage_red.txt', string_redpercentage,'\t');
    figure(1314)                %1314 is just the handle of the percentage_of_red drawing, it is just for avoid confusion
    plot(string_redpercentage,'rx-');
    xlabel('frame');
    ylabel('pertentage of red particles over all the particles');
    print(1314,'-dpng','percentage_red.png')
    close 1314
    cd('..');
end


%system('convert *.png -compose darken -flatten -quality 90 sum.jpg');

