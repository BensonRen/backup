% called by: trackcalall_Isobe_Main_BSseperate_Ben'
function []=trackcalall_Isobe_BSseperate_Ben(m)
%nseg=20;
nseg=10;

m=unwrap(m)
%msdcal(m); savetxtfig('msd', 'msd');
%histdx(m,'dim',[1 2]); 
%savetxtfig('histdx', 'histdx'); xlim([0 1]); screen2png('histdxA'); xlim([0 0.2]);   %screen2png('histdxB'); 
%system('gp output histdxsqr.png using "((column(1))**2):2" s histdx.txt l  xlabel "dx^2" ylabel "P(dx)"'); % to compare with normal distribution
%showdispall(m, 90000, '06_nseg2', 'nseg',2, 'drmin',0.6, 'dr2',10, 'MS',8);
%showdispall(m, 90000, '06', 'nseg',nseg, 'drmin',0.6, 'dr2',10, 'MS',8);
%showdispall(m, 90000, '0', 'nseg',nseg, 'drmin',0, 'dr2',10, 'MS',8);
%showdispall(m, 90000, '0_nseg500', 'nseg',500, 'drmin',0, 'dr2',10, 'MS',8);
%system('mv -f *-90000/* .; rmdir *-90000');
showdispall(m, 1000, '0', 'nseg',nseg, 'drmin',0, 'dr2',10, 'MS',8);
%showdispall_Ben(m, 1000, '0_show_radius','nseg',nseg, 'drmin',0,'dr2',10, 'MS',8,'show_radius');
showdispall(m, 300, '0_show_radius','nseg',nseg, 'drmin',0, 'dr2',10, 'MS',8,'show_radius');
showdispall(m, 100,'0_show_radius', 'nseg',nseg, 'drmin',0, 'dr2',10, 'MS',8,'show_radius');

%showdispall_Ben(m, 50, '0', 'nseg',nseg, 'drmin',0, 'dr2',10, 'MS',8);
%showdispall_Ben(m, 20, '0', 'nseg',nseg, 'drmin',0, 'dr2',10, 'MS',8);

w=findhopret(m,'Dt',1000,'dt',5)
