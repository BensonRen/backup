%%This program deal with large data sets from Isobe data, aligned by
%%particle. The input variables are:
%n= number of particles illustrated by the parameter.txt
%nframe= number of frames in total illustrated by the parameter.txt
%fget= number of frames you want to get in your o
filename='N0064_0720_AVE'
format='.dat'
file=[filename format];
data = dlmread(file);
[h, l] = size(data);
n=4096          %number of particles 
nframe=1000   %number of frames in total
fget=100     %the number of frames you want to extract
j=1;
datac=zeros(n*fget,l);
for i=1:nframe:h
    datac(j:j+fget-1,:)=data(i:i+fget-1,:);
    j=j+fget;
end
outname=[filename num2str(fget) format]
dlmwrite(outname, datac,'\t');