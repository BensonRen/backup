%%This function is to divide the Isobe data into big and small .dat files
%%to analyze them seperately
function [] = divide_BS_Isobe(filename)
file_0=[filename num2str(0)];       %construct a file_0 name for the dat0 
data=dlmread(filename);             %read the posisiton file
data0=dlmread(file_0);              %read the radius file
[h, l] = size(data);                %size of data
[h0,l0]= size(data0);               %size of data0
radius_small=min(data0(:,3));       %find the big and small radius of the particle
radius_big=max(data0(:,3));
big_identifier=(data0(:,3)==radius_big); %get the logical array of big particles
nbig=sum(big_identifier);
nsmall=sum(~big_identifier);
number_big_small=[nbig nsmall];
dlmwrite('number_big_small.txt',number_big_small,'\t');

n=h0          %number of particles 
nframe=h/h0   %number of frames in total
data_big=zeros(nbig*nframe,2);
data_small=zeros(nsmall*nframe,2);
j=1;
k=1;
for i=1:n:h
    temp=data(i:i+n-1,:);
    data_big(j:j+nbig-1,:)=temp(big_identifier,:);
    data_small(k:k+nsmall-1,:)=temp(~big_identifier,:);
    j=j+nbig;
    k=k+nsmall;
end
bigfile=[filename 'big'];
smallfile=[filename 'small'];
dlmwrite(bigfile,data_big,'\t');
dlmwrite(smallfile,data_small,'\t');
    
