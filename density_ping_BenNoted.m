clear; %clear workspace

%%Setting parameters
%define radius range
Rmin = 9;
Rmax = 30;

%define particle related parameters: Magnification and radius_ratio
magnification = 1.11
radius_ratio = 1.22;

%define the boundary of the 'within' region
a = 400
b = 1400
c = 200
d = 800

fig1 = figure();
%Open a new figure1

A = imread('../0016.jpg');
% A = imread('Image/00000.jpg');
%read the image according to the route given

imshow(A)
%display the image chosen

[centers, radii] = imfindcircles(A,[Rmin Rmax],'ObjectPolarity','bright');
%find all the bright circles in the image using built-in function
%imfindcircles and assign the centers and radius to the array [centers,
%radii]

ignore = (((centers(:,1)-magnification*radii)<=0) |  ((centers(:,1)+magnification*radii)>=1980) |  ((centers(:,2)-magnification*radii)<=0) | ((centers(:,2)+magnification*radii)>=1114));
%if the circle is out of the range of the size of the picture, identify it
%to be variable ignore

%delete those out-of-range circles
centers = centers(~ignore,:);
radii = radii(~ignore);

%assign x and y value of the particles
x = centers(:,1);
y = centers(:,2);

%make a area that is called 'within' which is designated by the chosen
%value at the start of the program, with boundaries abcd
within = (x>=a) & (x<b) & (y>=c) & (y<d);

%open figure2
fig2 = figure();

%make the histogram of the radius of all the particles inside the 'within'
%range chosen
h = histogram(radii(within),100); 
hold on

%%fit the shape of the histogram to a gaussian function
%take the middle value of each and every BinEdges
hist_bins = transpose((h.BinEdges(1:end-1)+h.BinEdges(2:end))./2);
%fit the mid_point and corresponding value with a two-term gaussian model
f = fit(hist_bins,transpose(h.Values),'gauss2');
plot(f);

%result f is a cfit-type object with parameters assigned:
%a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2)

%determine which is the big radius(center of the gaussians)
if f.b1 > f.b2
    radius_big = f.b1;
    radius_small = f.b2;
else
    radius_big = f.b2;
    radius_small = f.b1;
end;

%divide radius from small to big using 0.01 as interval
xg = (radius_small:0.01:radius_big);
yg = f.a1.*exp(-((xg-f.b1)./f.c1).^2) + f.a2.*exp(-((xg-f.b2)./f.c2).^2);
%find the fitted value 

%=====THIS CAN BE DENOTED SIMPLER USING THE FOLLOWING CODE:====
%yg=f(xg)';

% plot(xg,yg,'linewidth',3)

%Set the boudary as the x which is the lowest point 
boundary = xg(find(yg == min(yg)))
%plot the boundary line from min to max at the boundary point using thick
%line
plot([boundary, boundary], [min(yg), max(yg)],'linewidth',3);

%compute the radius ratio
original_radius_ratio = radius_big/radius_small
% radius_add = (radius_big - radius_small*radius_ratio)/(radius_ratio - 1);

%create a logical sequence determining which particle is a big particle
big = (radii>boundary);
%select those particles that is within the range 'within'
radii_within = radii(within);

%set var mole_ratio equal to number of particles that is larger that
%boundary to those small than boundary inside the region 'within'
mole_ratio = length(radii_within(radii_within>boundary)) / length(radii_within(radii_within<boundary))

%%% radii(~big) = (radii(~big) + radius_add) * radius_small / (radius_small+radius_add);
%%% radii(big) = (radii(big) + radius_add) * radius_small / (radius_small+radius_add);

figure(fig2);
h_after = histogram(radii(within),100);
legend('before','curve fitting','boundary','scaled');

%magnify radii of both big and small particles
radii(~big) = radii(~big) * magnification;
radii(big) = radii(big) * magnification;

%=====THE ABOVE TWO LINE CODE IS EQUIVALENT TO:====
%radii=radii*magnification;

%Plot all the circles using built-in function viscircles
figure(fig1);
%viscircles(centers(big,:), radii(big),'Color','b');
%viscircles(centers(~big,:), radii(~big),'Color','r');
viscircles(centers, radii);
hold on
%plot the boundaries of the 'within'
plot([a,b,b,a,a], [c,c,d,d,c], 'k-', 'linewidth', 3)

%Open a figure 3
fig3 = figure();
%use built-in function set to reverse the Y-axis, gca is the current axis
%handle, this is to make the circle plot the same as the original image
%since the image counts the y axis from the top to the bottom
set(gca,'Ydir','reverse')
%set x-y to be the same length
axis('equal')

%viscircles(centers(big,:), radii(big),'Color','b');
%viscircles(centers(~big,:), radii(~big),'Color','r');

%Draw all the particles
viscircles(centers, radii);
hold on

%color the centers of circle seperately
plot(centers(big,1), centers(big,2), 'b.');
plot(centers(~big,1), centers(~big,2), 'r.');
hold on
%draw the bounday
plot([a,b,b,a,a], [c,c,d,d,c], 'k-', 'linewidth', 3)

%open figure4 to draw the particles that is in side the 'within' area
fig4 = figure();
set(gca,'Ydir','reverse')
axis('equal')
%viscircles(centers(within,:), radii(within),'Color','b');
%viscircles(centers(~within,:), radii(~within),'Color','r');
viscircles(centers, radii);
hold on
plot(centers(within,1), centers(within,2), 'b.');
plot(centers(~within,1), centers(~within,2), 'r.');
hold on
plot([a,b,b,a,a], [c,c,d,d,c], 'k-', 'linewidth', 3)


%%compute the density of the 'within' region
s = sum(radii(within) .^ 2 *pi);
density = s/(b-a)/(d-c)

%%save the figures generated
saveas(fig1,'overlay.fig')
saveas(fig2,'distribution.fig')
saveas(fig3,'particles.fig')
saveas(fig4,'particles_within.fig')

saveas(fig1,'overlay.png')
saveas(fig2,'distribution.png')
saveas(fig3,'particles.png')
saveas(fig4,'particles_within.png')

T1 = table;
T1.Parameter = {'Radius ratio'; 'Average small radius'; 'a'; 'b'; 'c'; 'd'; 'magnification'};
T1.Value = [radius_ratio; radius_small; a; b; c; d; magnification];
writetable(T1, 'parameter_radius.txt', 'Delimiter', '\t');

T2 = table;
T2.Parameter = {'Radius ratio'; 'Total number of particles'; 'Mole fraction (big to small)'; 'Original radius ratio'; 'magnification'; 'Density'};
T2.Value = [radius_ratio; length(radii); mole_ratio; original_radius_ratio; magnification; density];
writetable(T2, 'parameter_density.txt', 'Delimiter', '\t');
