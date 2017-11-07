
drawnow update

n = length( ref )/3;
eqpos = reshape(ref, 3, n)';
eqpos = [x0.'; eqpos(2:M+1,:); eqpos(1,:)];

figure(1), set(gcf, 'Color','white');
clf

subplot(3,3,[1:6]);
TOL = 0.00;
p = patch([-0.2, 1.2, 1.2, -0.2], [WALL-TOL, WALL-TOL, WALL-TOL, WALL-TOL], [-4, -4, 1, 1], 'g');
hold on;

plot3(eqpos(:, 1), eqpos(:, 2), eqpos(:, 3), '-ob', ...
    'MarkerSize', 7.5, 'MarkerFaceColor', 'b', 'linewidth', 0.2);

view([-135 45*3/4]);

xlim([-0.2 1.2]);
ylim([-0.2 1.2]);
zlim([-4 1]);

grid on;

set(gca, 'Box', 'on');

curpos = reshape(state_sim(end,:), 3, n)';
curpos = [x0.'; curpos(2:M+1,:); curpos(1,:)];

plot3(curpos(:, 1), curpos(:, 2), curpos(:, 3), '-or', ...
    'MarkerSize', 7.5, 'MarkerFaceColor', 'r', 'linewidth', 0.2);

xlabel( 'x [m]' );
ylabel( 'y [m]' );
zlabel( 'z [m]' );


plot3(eqpos(1, 1), eqpos(1, 2), eqpos(1, 3), '-ok', ...
    'MarkerSize', 7.5, 'MarkerFaceColor', 'k', 'linewidth', 0.2);


if( ~isempty(controls_MPC) )
    subplot(3,3,7);
    stairs(time(1:end-1),controls_MPC(:,1));
    grid on;
    xlabel('t [s]','FontSize',13);    ylabel('u1','FontSize',13)
    
    subplot(3,3,8);
    stairs(time(1:end-1),controls_MPC(:,2));
    grid on;
    xlabel('t [s]','FontSize',13);    ylabel('u2','FontSize',13)
    
    subplot(3,3,9);
    stairs(time(1:end-1),controls_MPC(:,3));
    grid on;
    xlabel('t [s]','FontSize',13);    ylabel('u3','FontSize',13)
end
