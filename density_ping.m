clear;                          %clear the workspace
Rmin = 9;                      %Set some parameters
Rmax = 20;
sen = 0.88;
bandpass = 1;
% threshold = 0.8;
% magnification = 1.11
small_add = 3.1
big_add = 3.1
radius_add = small_add
% radius_ratio = 1.22;

fig1 = figure();                %open new figure
frame = '0152902.jpg';           %set the frame to be analyze
A = imread(frame);              %read frame
% imagesc(A)          
A = rgb2gray(A);                %change to grad scale
A = bpass(A, bandpass, 40);     %use bpass to filter
% A = imgaussfilt(A, 2);
% A = im2bw(A, threshold);
colormap('gray');               %set the color map of current figure 
imagesc(A)                      %show image
% imshow(A)

axis('equal')

% pk = pkfnd(A, 25, 29);
% cnt = cntrd(A, pk, 29);
% metric = cnt(:,3);
% centers = cnt(:, 1:2);
% radii = cnt(:,4).^0.5;
% 
% metric_threshold = 5000;
% centers = centers(metric>metric_threshold,:);
% radii = radii(metric>metric_threshold);
% metric = metric(metric>metric_threshold);

[centers, radii, metric] = imfindcircles(A,[Rmin Rmax],'ObjectPolarity','bright','Sensitivity',sen);
%assign the x,y value 
x = centers(:,1);
y = centers(:,2);
a = 400
b = 1500
c = 100
d = 900
% [l, w] = size(A);
% a = w*0.1; b = w*0.9; c = l*0.15 ;d = l*0.85;

%make a area that is called 'within' which is designated by the chosen
%value at the start of the program, with boundaries abcd
within = (x>=a) & (x<b) & (y>=c) & (y<d);

% fig5 = figure();
% h = histogram(metric(within),100);

%open figure2
fig2 = figure();

%make the histogram of the radius of all the particles inside the 'within'
%range chosen
% h = histogram(radii,100); 
h = histogram(radii(within), 100); 
hold on

%%fit the shape of the histogram to a gaussian function
%take the middle value of each and every BinEdges
hist_bins = transpose((h.BinEdges(1:end-1)+h.BinEdges(2:end))./2);

%fit the mid_point and corresponding value with a two-term gaussian model
f = fit(hist_bins,transpose(h.Values),'gauss2');
plot(f);

%determine which is the big radius(center of the gaussians)
if f.b1 > f.b2
    radius_big = f.b1
    radius_small = f.b2
else
    radius_big = f.b2
    radius_small = f.b1
end;

%divide radius from small to big using 0.01 as interval
xg = (radius_small:0.01:radius_big);
% yg = f.a1.*exp(-((xg-f.b1)./f.c1).^2) + f.a2.*exp(-((xg-f.b2)./f.c2).^2);
yg = f(xg)';
% plot(xg,yg,'linewidth',3)

%Set the boudary as the x which is the lowest point 
boundary = xg(find(yg == min(yg)))
%plot the boundary line from min to max at the boundary point using thick
%line
plot([boundary, boundary], [min(yg), max(yg)],'linewidth',3);

%compute the radius ratio
original_radius_ratio = radius_big/radius_small
radius_ratio = (radius_big + big_add)/(radius_small + small_add)
add_1_4 = (radius_big-1.4*radius_small)/0.4

% radius_add = (radius_big - radius_small*radius_ratio)/(radius_ratio - 1);

%create a logical sequence determining which particle is a big particle
big = (radii>boundary);
%select those particles that is within the range 'within'
radii_within = radii(within);

%Set the particle numbers
small_particle_number = length(radii_within(radii_within<boundary))
big_particle_number = length(radii_within(radii_within>boundary))
big_particle_number_ratio = big_particle_number / length(radii_within)

% radii(~big) = (radii(~big) + radius_add) * radius_small / (radius_small+radius_add);
% radii(big) = (radii(big) + radius_add) * radius_small / (radius_small+radius_add);
radii(~big) = radii(~big) + small_add;
radii(big) = radii(big) + big_add;

% figure(fig2);
% % h_after = histogram(radii,100);
% h_after = histogram(radii(within), 100);
% legend('before','curve fitting','boundary','scaled');

% radii(~big) = radii(~big) * magnification;
% radii(big) = radii(big) * magnification;

%set the big and small particles
radii_big = radii(big);
radii_small = radii(~big);
centers_big = centers(big,:);
centers_small = centers(~big,:);

figure(fig1);
% for i = (1:length(radii_big))
%     r = radii_big(i);
%     cen = centers_big(i,:);
%     pos = [cen-r 2*r 2*r];
%     rec = rectangle('Position', pos, 'Curvature', [1 1], 'FaceColor', 'b', 'Edgecolor', 'none');
%     alpha(rec, 0.1)
% end;

%draw the big particles with blue
viscircles(centers_big, radii_big,'EdgeColor','b', 'linewidth', 1);
% for i = (1:length(radii_small))
%     r = radii_small(i);
%     cen = centers_small(i,:);
%     pos = [cen-r 2*r 2*r];
%     rec = rectangle('Position', pos, 'Curvature', [1 1], 'FaceColor', 'r', 'Edgecolor', 'none');
% end;

%draw the small particles with red
viscircles(centers_small, radii_small,'EdgeColor','r', 'linewidth', 1);
hold on
%draw the within area
plot([a,b,b,a,a], [c,c,d,d,c], 'k-', 'linewidth', 3)
hold on
% for i = (1:length(radii))
%     text(centers(i,1), centers(i,2), num2str(i));
%     hold on
% end

fig3 = figure();
set(gca,'Ydir','reverse')
axis('equal')
% for i = (1:length(radii_big))
%     r = radii_big(i);
%     cen = centers_big(i,:);
%     pos = [cen-r 2*r 2*r];
%     rectangle('Position', pos, 'Curvature', [1 1], 'FaceColor', 'b', 'Edgecolor', 'none')
% end;
viscircles(centers_big, radii_big,'EdgeColor','b', 'linewidth', 1);
% for i = (1:length(radii_small))
%     r = radii_small(i);
%     cen = centers_small(i,:);
%     pos = [cen-r 2*r 2*r];
%     rectangle('Position', pos, 'Curvature', [1 1], 'FaceColor', 'r', 'Edgecolor', 'none')
% end;
% hold on;
viscircles(centers_small, radii_small,'EdgeColor','r', 'linewidth', 1);
hold on;
plot(centers_big(:,1), centers_big(:,2), 'b.');
plot(centers_small(:,1), centers_small(:,2), 'r.');
hold on
plot([a,b,b,a,a], [c,c,d,d,c], 'k-', 'linewidth', 3)


% fig4 = figure();
% set(gca,'Ydir','reverse')
% axis('equal')
% viscircles(centers(within,:), radii(within),'Color','b', 'linewidth', 1);
% viscircles(centers(~within,:), radii(~within),'Color','r', 'linewidth', 1);
% hold on
% plot(centers(within,1), centers(within,2), 'b.');
% plot(centers(~within,1), centers(~within,2), 'r.');
% hold on
% plot([a,b,b,a,a], [c,c,d,d,c], 'k-', 'linewidth', 3)

%compute the density within the region
s = sum(radii(within) .^ 2 *pi);
density = s/(b-a)/(d-c)


% saveas(fig1,'overlay.fig')
% saveas(fig2,'distribution.fig')
% saveas(fig3,'particles.fig')
% saveas(fig4,'particles_within.fig')
% 
% saveas(fig1,'overlay.png')
% saveas(fig2,'distribution.png')
% saveas(fig3,'particles.png')
% saveas(fig4,'particles_within.png')
% 
T = table;
T.Parameter = {'Density'; 'Average small radius'; 'a'; 'b'; 'c'; 'd'; 'Radius ratio'; 'Original radius ratio'; 'Total number of particles'; 'Big particle fraction (big to total)'; 'Frame'; 'Small_radius_add'; 'Big_radius_add'; 'Radius_min'; 'Radius_max'; 'Sensitivity'; 'Bandpass'};
T.Value = {density; radius_small+radius_add; a; b; c; d; radius_ratio; original_radius_ratio; length(radii); big_particle_number_ratio; frame; small_add; big_add; Rmin; Rmax; sen; bandpass};
writetable(T, 'parameter_radius.txt', 'Delimiter', '\t');
