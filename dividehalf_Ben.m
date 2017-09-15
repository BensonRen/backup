%%DIVIDEHALLF_BEN is a script that takes in a 0.dat file and take every
%%other frame of the particle to divide the 0.dat into half of the original
%%size, of course, lose of information is unavoidable but it gives up
%%capability to handle super large data sets and see the long time frame
%%picture

data_original=dlmread('0.dat'); %read the data in 0.dat
fre_table=tabulate(data_original(:,6));     %tabulate the frequency of partilces
np=max(data_original(:,6));     %number of particles
nf=max(data_original(:,5));     %number of frames
ncol=size(data_original,1);     %number of coloums of data
datac=zeros(fix(ncol/2)+np,6);  %pre-allocation of datac
j=1;k=1;
%%assign every other frame into each particle
for i=1:np  %for each particle
    nfi=fre_table(i,2);
    if mod(nfi,2)==0     %number of frame is even
        datac(k:k+nfi/2-1,:)=data_original(j:2:j+nfi-1,:);
        k=k+nfi/2;
    else    %numer of frames is odd
        datac(k:k+(nfi-1)/2,:)=data_original(j:2:j+nfi-1,:);
        k=k+(nfi+1)/2;
    end
    j=j+nfi;     %adder to j
end

%%eliminate the extra 0s
for i=fix(ncol/2)-np:fix(ncol/2)+np
    if isequal(datac(i,1:4),[0 0 0 0])
        check=i-1;
        break
    end
end
datac(:,5)=round(datac(:,5)/2);
datac=datac(1:check,:);
dlmwrite('c.dat', datac,'\t');    %save the file into c.dat




