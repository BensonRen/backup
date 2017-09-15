%ADDFRAMEALL takes some .dat file and combine them together in the sequence
%of their prefix, e.g. 1.dat+2.dat+3.dat etc
%It calls another function called addframe_Ben to finish the concatnation
%of the 0.dats
%One have to manually put the 0.dats together inside the folder to work
%with it
list = dir;              %Get all the files in this
n=0;

%Check how many .dat files there are to determine the number of times to do
%the addframe
for i=1:length(dir)                 %for every file in current directory
    k=strfind(list(i).name,'.dat'); %find '.dat' inside the file name
    if (~isempty(k))                %if k is not empty, which means it is .dat file
        n=n+1;                      %adder 1 to the #dat files
    end
end

for i=3:n                           %from 3.dat
    addframe_Ben;                   %addframe
    status=unix('mv c.dat 1.dat');  %change c.dat to 1.dat
    file2=[num2str(i) '.dat'];      %make the name of i.dat
    status=unix(['mv ' file2 ' 2.dat']); %change it to 2.dat
end
addframe_Ben;
status=unix('mv c.dat 0.dat');      %make the final result to 0.dat
status=unix('cp -r ../density .');  %get the outside dentisy folder in

