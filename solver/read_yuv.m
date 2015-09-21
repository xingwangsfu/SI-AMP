
function img = read_yuv(path, row, col)
fid1 = fopen(path,'r');

% row = 1024;
% col = 768;

img= fread(fid1, [row,col],'uchar');

% fid2 = fopen('C:\Users\xingw\Desktop\wx\output_balloons\balloons_03_synth.yuv.all.yuv','r');
% 
% T_inter =  fread(fid2, [row,col],'uchar');