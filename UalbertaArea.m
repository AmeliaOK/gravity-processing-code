function[] = UalbertaArea(Transparency)
%% GIVES THE IMAGE OF DOROTHY PROJECTED OVER THE DATA BUT TRANSPARENT [0:1]
hold on;
A = imread('UofAmap.jpg');
B = flip(A);
% R = imref2d(size(B),[332185 333374],[5933521 5934288]); 
%R = imref2d(size(B),[332179 333368],[5933509 5934276]); 

X1= 332362.635555076;
X2= 332863.687513639;
Y1= 5934257.56204513;
Y2= 5933605.0394435;

R = imref2d(size(B),[X1 X2],[Y2 Y1]); 


axis on;
im = imshow(B,R);
im.AlphaData = Transparency;
set(gca, 'YDir', 'normal');
hold on;
end
