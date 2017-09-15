%function [] = tracking_ping(Radius_add, rmin, rmax)
function [] = tracking_ping(Rmin, Rmax)
%function [] = tracking_ping(Rmin, Rmax, threshold)
warning('off','all');

% magnification = 1.11
%radius_add = 2
%Rmin = 8
%Rmax = 30
%radius_add = Radius_add
Rmin
Rmax
%threshold

%read all the parameters into the table par
par = readtable('../density/parameter_radius.txt','Delimiter','\t');
radius_big_real = par.Value(2) * par.Value(7)

%detect if it exist a datafile (Or is it continuing a previous job)
if ~exist('0000000.dat') %This is a new job, start from the begining
    datafile = [];
    %use built-in function dlmwrite to write matrix to ASCII-delimited file
    dlmwrite('0000000.dat', [], '\t');
    first_frame = 1;
else %This is not a new job, start from the end of the datafile and keep doing
    datafile = dlmread('0000000.dat');
    first_frame = datafile(end, end) + 2;
end

%set up some directory names to access and determine the numebr of frames
frames_dir = 'frame/'
frames = dir ( fullfile(frames_dir,'*.jpg') );
frames_num = length(frames)

%fig = figure(1);

for f = (first_frame:frames_num)
    %for each frame, imread it and convert it into gray scale
    A = imread(fullfile(frames_dir,frames(f).name));
    A = rgb2gray(A);
    
    %%%%=======What do this bpass do exactly?
    A = bpass(A, 2, radius_big_real);
    %A = im2bw(A, 0.8);
    %     imshow(A)
    %=======================================
    
    %read the parameters of the particles using built-in function
    %imfindcircles
    [centers, radii] = imfindcircles(A,[Rmin Rmax],'ObjectPolarity','bright','Sensitivity',0.87);
    
    %h = viscircles(centers, radii);
    %         h = circle(centers(:,1), centers(:,2), radii, 'b');
    axis('equal');
    drawnow;
    %         saveas(fig,['particles', num2str(f), '.png'])
    %frame = getframe(1);
    %im{f} = frame2im(frame);
    
    
    %%%%==============ADD THIS WOULD MAKE BETTER PERFORMANCE============
    %temp=zeros(length(centers),5);
    %temp(:,1:2)=centers;
    %===================================================================
    
    %use a temp array to store the position of the circles
    temp = centers;
    
    
    
    %mass = pi .* ((radii + radius_add) .^ 2);
    
    %compute the area of each particle assign it to mass, then assign it to
    %the third colomn of temp
    mass = pi .* (radii .^ 2);
    temp(:,3) = mass;
    
    %     temp(:,4) = radii .* magnification;
    %temp(:,4) = radii + radius_add;
    
    %write the radius of the particles as 4th colomn of temp
    temp(:,4) = radii;
    %write the frame number to the 5th colomn of temp
    temp(:,5) = f-1;
    %write temp into the file called '0000000.dat'
    dlmwrite('0000000.dat',temp,'-append','delimiter','\t');
    %append temp to the end of the datafile
    datafile = [datafile; temp];
    %print our the number of particles found in this frame
    disp([fullfile(frames_dir, frames(f).name) ': ' int2str(length(centers)) ' particles found']);
    
    %delete(h);
    %         for i = 1:length(h)
    %             delete(h{i})
    %         end
end

%gif_name = 'test.gif';
%for idx = 1:frames_num
%    [A, map] = rgb2ind(im{idx}, 256);
%    if idx == 1
%        imwrite(A, map, gif_name, 'gif', 'LoopCount',Inf,'DelayTime',1);
%    else
%        imwrite(A ,map, gif_name,'gif','WriteMode','append','DelayTime',1);
%    end
%end

%     fig2 = figure();
%radii_allframe = datafile(:,4);
%hist = histogram(radii_allframe,100);
%     hold on
%hist_bins = transpose((hist.BinEdges(1:end-1)+hist.BinEdges(2:end))./2);
%f = fit(hist_bins,transpose(hist.Values),'gauss2');
%     plot(f);

%if f.b1 > f.b2
%    radius_big = f.b1
%    radius_small = f.b2
%else
%    radius_big = f.b2
%    radius_small = f.b1
%end;

%Setting the parameters for the function track.m 
param.mem = 1; %number of time steps that is allowed for a particle to be missing
param.good = 0.2 * frames_num %%============Can't understand this==========
param.dim = 2 %%============Can't understand this==========
param.quiet = 1 %no text output

%data_tracked = track(datafile ,radius_big, param);
%data_tracked = track(datafile ,max(datafile(:,4)), param);
data_tracked = track(datafile ,radius_big_real, param);
%data_tracked = track(datafile ,20 , param);
data_tracked(:,end) = data_tracked(:,end) - 1;
dlmwrite('0.dat', data_tracked,'\t');
%delete(hist);
%unix('nohup matlab -r trackcalall &');
%exit;
