% Initialize ballgun on the LEFT = 0, CENTER = 1, RIGHT = 2
% with random ball positions and velocities
% Gaussian r.v. with standard deviation std

function [ball_pos,ball_vel] = initBallGun(std,side)

    loadTennisTableValues();
	ballgun = [table_center; ...
               dist_to_table - table_length - 0.2; ...
               floor_level + table_height + 0.15];
	switch side 
        case 'LEFT'
            disp('Setting ballgun to left side');
            good_ball_vel = [-1.08; 4.00; 3.84];
            ballgun(1) = ballgun(1) + 0.4;
        case 'CENTRE'
            disp('Setting ballgun to center');
            good_ball_vel = [0.0; 4.0; 3.84];
        case 'RIGHT'
            disp('Setting ballgun to right side');
            good_ball_vel = [1.08; 4.80; 3.84];
            ballgun(1) = ballgun(1) - 0.4;
        otherwise
            error('Ballgun location not identified!');
            % do nothing ballgun is already centred.
    end
	rand_ball_pos = ballgun + std * randn(3,1);
	rand_ball_vel = good_ball_vel + std * rand(3,1);

	ball_pos = rand_ball_pos;
	ball_vel = rand_ball_vel;


end

