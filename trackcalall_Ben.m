function trackcalall_Ben(m)
nseg=20;
%m=unwrap(m);
%w=findhopret(m)
msdcal(m); savetxtfig('msd', 'msd');
histdx(m,'dim',[1 2]); 
savetxtfig('histdx', 'histdx'); xlim([0 1]); screen2png('histdxA'); xlim([0 0.2]); screen2png('histdxB'); 
system('gp output histdxsqr.png using "((column(1))**2):2" s histdx.txt l  xlabel "dx^2" ylabel "P(dx)"'); % to compare with normal distribution
showdispall_Ben(m, 90000, '06_nseg2', 'nseg',2, 'drmin',0.6, 'dr2',10, 'MS',8);
showdispall_Ben(m, 90000, '06', 'nseg',nseg, 'drmin',0.6, 'dr2',10, 'MS',8);
showdispall_Ben(m, 90000, '0', 'nseg',nseg, 'drmin',0, 'dr2',10, 'MS',8);
%showdispall_Ben(m, 90000, '0_nseg500', 'nseg',500, 'drmin',0, 'dr2',10, 'MS',8);
system('mv -f *-90000/* .; rmdir *-90000');
%showdispall_Ben(m, 10000, '0', 'nseg',nseg, 'dr2',10, 'MS',8);
%showdispall_Ben(m, 10000, '0_show_radius','nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
%showdispall_Ben(m, 10000, '0_show_radius_real_cen','nseg',20, 'dr2',10, 'MS',8,'show_radius','drpath','real_cen');

%showdispall_Ben(m, 10000, '0_top_radius', 'top_radius', 10,'nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
showdispall_Ben(m, 5000, '0', 'nseg',nseg, 'dr2',10, 'MS',8);
%showdispall_Ben(m, 5000, '0_show_radius','nseg',nseg, 'dr2',10, 'MS',8,'show_radius','gif');
showdispall_Ben(m, 5000, '0_show_radius','nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
%showdispall_Ben(m, 5000, '0_show_radius_real_cen','nseg',20, 'dr2',10, 'MS',8,'show_radius','drpath','real_cen');

%showdispall_Ben(m, 5000, '0_top_radius', 'top_radius', 10,'nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
showdispall_Ben(m, 2000, '0', 'nseg',nseg, 'dr2',10, 'MS',8);
%showdispall_Ben(m, 2000, '0_show_radius','nseg',nseg, 'dr2',10, 'MS',8,'show_radius','gif');
showdispall_Ben(m, 2000, '0_show_radius','nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
%showdispall_Ben(m, 2000, '0_show_radius_real_cen','nseg',nseg, 'dr2',10, 'MS',8,'show_radius','drpath','accurate','real_cen');
%showdispall_Ben(m, 2000, '0_top_radius', 'top_radius', 10,'nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
showdispall_Ben(m, 1000, '0','nseg',nseg, 'dr2',10, 'MS',8);
showdispall_Ben(m, 1000, '0_show_radius','nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
%showdispall_Ben(m, 1000, '0_show_radius_real_cen','nseg',nseg, 'dr2',10, 'MS',8,'show_radius','drpath','accurate','real_cen');
%showdispall_Ben(m, 1000, '0_top_radius', 'top_radius', 10,'nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
%showdispall_Ben(m, 300, '0','nseg',nseg , 'dr2',10, 'MS',8);    
%showdispall_Ben(m, 300, '0_show_radius','nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
%showdispall_Ben(m, 300, '0_show_radius_real_cen','nseg',nseg, 'dr2',10, 'MS',8,'show_radius','drpath','accurate','real_cen');
%showdispall_Ben(m, 300, '0_top_radius', 'top_radius', 10,'nseg',nseg, 'dr2',10, 'MS',8,'show_radius');

%showdispall_Ben(m, 100, '0','nseg',nseg, 'dr2',10, 'MS',8);
%showdispall_Ben(m, 100, '0_show_radius','nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
%showdispall_Ben(m, 100, '0_show_radius_real_cen','nseg',5, 'dr2',10, 'MS',8,'show_radius','drpath','accurate','real_cen');
%showdispall_Ben(m, 100, '0_top_radius', 'top_radius', 10, 'nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
%showdispall_Ben(m, 30, '0','nseg',5, 'dr2',10, 'MS',8);
%showdispall_Ben(m, 30, '0_show_radius','nseg',nseg, 'dr2',10, 'MS',8,'show_radius', 'drpath','accurate');
%showdispall_Ben(m, 30, '0_show_radius_real_cen','nseg',5, 'dr2',10, 'MS',8,'show_radius', 'drpath','accurate','real_cen');
%showdispall_Ben(m, 30, '0_top_radius', 'top_radius', 10, 'nseg',nseg, 'dr2',10, 'MS',8,'show_radius');
%showdispall_Ben(m, 10, '0','nseg',5, 'dr2',10, 'MS',8);
%showdispall_Ben(m, 10, '0_show_radius','nseg',nseg, 'dr2',10, 'MS',8,'show_radius', 'drpath','accurate');
%showdispall_Ben(m, 10, '0_show_radius_real_cen','nseg',5, 'dr2',10, 'MS',8,'show_radius', 'drpath','accurate','real_cen');

w=findhopret(m, 'Dt', 5000)

%%TRACKCALALL FOR BIG PARTICLES
%exit;

end