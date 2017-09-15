function [] = tracking_ping_all()
%function [] = tracking_ping_all(rmin, rmax)
%function [] = tracking_ping_all(Radius_add, rmin, rmax)

%set all the warning states to be off ?Why
warning('off','all'); 

run_max = 3

% magnification = 1.11
% radius_add = 2
%Rmin = 8
%Rmax = 30
%radius_add = Radius_add
%Rmin = rmin
%Rmax = rmax

%read the parameter_radius file as a table using \t as seperater
par = readtable('density/parameter_radius.txt','Delimiter','\t');
Rmin = par.Value(end-1)
Rmax = par.Value(end)
threshold = 0.8

%make the current directory as a struct, which contains the folder
%names,date,bytes, datenum(Date and time serial number) and a logical 
%indicator that signifies it is a dir or not
list = dir;
%Make the folder '.'(Current folder) and '..'(Upper folder) not a directory
for i = 1:length(list)
    if strcmp(list(i).name, '.') || strcmp(list(i).name, '..')
        list(i).isdir = false;
    end
end

%construct a struct that only contains the folders inside
folders = list([list.isdir]);

%define an empty cell called run_now_dir, with run_max number of cells
run_now_dir = cell(1,run_max)
run_now = 0

%run the folder one by one using i as counter
for i = length(folders):-1:1
    %when run_now==run_max, check if the 0.dat has been generated, is so,
    %run_now--, continue doing tracking
    while run_now == run_max
        pause(300)
        for j = 1:run_max
            if exist(fullfile(run_now_dir{j},'0.dat'))
                run_now_dir{j} = []
                run_now = run_now - 1
                break
            end
        end
    end
    
    %get the current folder name
    folder_name = folders(i).name
    frames_dir_name = 'frame';
    %use a built-in function fullfile to build the full path of the file
    frames_dir = fullfile(folder_name,frames_dir_name);
    
    %check if there is a 'frame' folder inside all the folders, if yes then
    %cd to that folder
    if ~exist(frames_dir) || exist(fullfile(folder_name, '0.dat'))
        continue
    else
        cd(folder_name)
    end
    
    %usd the built-in function UNIX to execute unix scripts, which in here
    %is to call tracking_ping function
    
    %unix(['nohup matlab -r ''tracking_ping(' num2str(radius_add) ',' num2str(Rmin) ',' num2str(Rmax) ')'''  '&'  '>'  'nohup.out&']);
    unix(['nohup matlab -r ''tracking_ping('  num2str(Rmin) ',' num2str(Rmax) ',' num2str(threshold) ')'''  '&'  '>'  'nohup.out&']);
    %unix(['nohup matlab -r ''tracking_ping('  num2str(Rmin) ',' num2str(Rmax) ')'''  '&'  '>'  'nohup.out&']);
    
    %assign the current folder name to run_now_dir{j}, run_now++
    %run_now denote the number of tracking that is currently doing
    %for safety of RAM, do not exceed run_max=3
    for j = 1:run_max
        if isempty(run_now_dir{j})
            run_now_dir{j} = folder_name
            run_now = run_now + 1
            break
        end
    end
    
    
    %return to the previous folder to prepare for the next iteration
    cd('..')
end
exit;
