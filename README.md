[![View Determine the distance between two ellipses (in 3D) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/26329-determine-the-distance-between-two-ellipses-in-3d)

[![Donate to Rody](https://i.stack.imgur.com/bneea.png)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4M7RMVNMKAXXQ&source=url)

# FEX-distanceEllipseEllipse

The problem of finding the geometric (minimum) distance between two arbitrary ellipses is surprisingly difficult. The general problem of finding all stationary points (minimum/maximum/saddle, no less than 12 possible points) has indeed been solved, but that algorithm is horrificly complex and would require thousands if not millions of operations once implemented.
The function distanceEllipseEllipse() is based on a somewhat more practical algorithm which limits itself to minima in the distance function. It is based on repeatedly finding the distance between a point and an ellipse, which can be done analytically (see my other post, distanceEllipsePoints.m). The algorithm by itself is not very robust (it frequently finds a local minimum, which is *not* the true distance).
This function executes the algorithm 4 times, for 4 different initial values, which greatly improves its robustness. Some numerical experimentation (comparing with a brute-force search) has shown that the true minimum distance is returned in more than 95% of the cases. The algorithm is implemented such that MATLAB's JIT-accelerator can accelerate it to the fullest extent, which makes it pretty fast and well-suited to handle large datasets requiring this calculation.
This is an implementation of the algorithm described in
Ik-Sung Kim: "An algorithm for finding the distance between two
% ellipses". Commun. Korean Math. Soc. 21 (2006), No.3, pp.559-567.

A copy-pastable example (also in header of M-file):

% Ellipse1 Ellipse2(=circle)
a = [2.0 1.0];
b = [0.5 1.0];
c = {[0,0,0], [-2,2,0]}; % location of centers
u = {[1,0,0], [1,0,0]}; % both oriented in XY-plane
v = {[0,1,0], [0,1,0]}; % to visualize them more easily

% plot the ellipses
f = 0:0.01:2*pi;
E1 = [a(1)*cos(f) + c{1}(1); b(1)*sin(f) + c{1}(2)];
E2 = [a(2)*cos(f) + c{2}(1); b(2)*sin(f) + c{2}(2)];
figure, hold on
plot(E1(1,:),E1(2,:),'r', E2(1,:),E2(2,:),'b')
axis equal

% run routine
[min_dist, fp_min, fs_min] = ...
distanceEllipseEllipse(a,b,c,u,v)

% plot the minimum distance returned
x = [a(1)*cos(fp_min) + c{1}(1), a(2)*cos(fs_min) + c{2}(1)];
y = [b(1)*sin(fp_min) + c{1}(2), b(2)*sin(fs_min) + c{2}(2)];
line(x,y,'color', 'k')

This should generate the screen shot given here.

If you like this work, please consider [a donation](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4M7RMVNMKAXXQ&source=url).
