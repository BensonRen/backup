
edi%%TRACKCALALL_BEN_MAIN is a m.script that takes the 0.dat and read into
%%three seperate struct:m, mbig, msmall and then track the trajectory of
%%the three
%it create three folders to store the image generated, all, big, small
%excecutor should note that this program needs a density folder in this
%directory and outside it for functions inside to call them both
nseg=20;
[m,mbig,msmall]=trackread_Ben_BSseperate('0.dat');
mkdir all
mkdir big
mkdir small
cd all
trackcalall_Ben(m);
m=[]
cd ../big
trackcalall_Ben(mbig);
mbig=[]
cd ../small
trackcalall_Ben(msmall);
msmall=[]
cd ..
quit