clear;clc;
filepath = 'C:\Users\A\Desktop\2023.02\STORM\20230314';

n = 3022;
for i = 1:n
    i
    I(:,:,i) = imread([filepath,'\a3_',num2str(i,'%.5d'),'.tif']);
end
options.color     = false;
options.compress  = 'no';
options.message   = true;
options.append    = false;
options.overwrite = true;
options.big       = true;
res = saveastiff(I, 'C:\Users\A\Desktop\1.tif', options)
