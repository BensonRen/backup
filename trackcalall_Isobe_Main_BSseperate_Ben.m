%%TRACKCALALL_Isobe_Mian_BSseperate is a m.script that takes two .dat file and read into
%%three seperate struct:m, mbig, msmall and then track the trajectory of
%%the three
%it create three folders to store the image generated, all, big, small
%excecutor should note that this program needs a density folder in this
%directory and outside it for functions inside to call them both
nseg=10
[m,mbig,msmall]=trackread_Isobe_BSseperate('N0064_0720_AVE.dat');
mkdir all
mkdir big
mkdir small
cd all
trackcalall_Isobe_BSseperate_Ben(m);
m=[]
cd ../big
trackcalall_Isobe_BSseperate_Ben(mbig);
mbig=[]
cd ../small
trackcalall_Isobe_BSseperate_Ben(msmall);
msmall=[]
cd ..
quit