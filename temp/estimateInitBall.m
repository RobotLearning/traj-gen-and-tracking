% Estimate the ball's initial position and velocity

function [ballInit,ballTime] = estimateInitBall(ballTime) 

	static  startflag= 1;
	static  threPos= 0.005; % m
	static  maxLen= 12;
	static  ballPos(12+1,4+1); % x,y,z,t
	static  temp1,temp2,num;
	static  numBallObservations =  0;
	static  PredictMode =  1;
	static  temp_time3, temp_time4;

	if (firsttime)
	
 	 	startflag = 0;
		temp_time3 = getTime();
	end
	temp_time4 = getTime();

	% dynamical saving window
	if (blobs(3).status == 1 && blobs(3).blob.x(2) > -2.8 && blobs(3).blob.x(2) < -1.2) || ...
	    (blobs(1).status == 1 && blobs(1).blob.x(2) > -2.0 && blobs(1).blob.x(2) < 0.5)
	
		if (blobs(3).status= = 1) 
            num = 3;
        else
            num = 1;
        end

		if (numBallObservations == 0)
		
			numBallObservations = numBallObservations+1;
			ballPos(numBallObservations,1) = blobs(num).blob.x(1); % m
			ballPos(numBallObservations,2) = blobs(num).blob.x(2); % m
			ballPos(numBallObservations,3) = blobs(num).blob.x(3); % m
			ballPos(numBallObservations,4) = (temp_time4-temp_time3)/1000/1000;% s
		else if(fabs(blobs(num).blob.x(2)-ballPos(numBallObservations,2)) <= threPos)
		
			%return;
		else if(fabs(blobs(num).blob.x(2)-ballPos(numBallObservations,2))>threPos)
		
			numBallObservations= numBallObservations+1;
			if (numBallObservations <= maxLen)
			
				ballPos(numBallObservations,1) = blobs(num).blob.x(1);
				ballPos(numBallObservations,2) = blobs(num).blob.x(2);
				ballPos(numBallObservations,3) = blobs(num).blob.x(3);
				ballPos(numBallObservations,4) = (temp_time4-temp_time3)/1000/1000;
			end

			if (numBallObservations>maxLen)
			
				for (temp1= 2;temp1 <numBallObservations;temp1= temp1+1)
				    for (temp2= 1;temp2 <= 4;temp2= temp2+1)
                        ballPos(temp1-1,temp2)= ballPos(temp1,temp2);
                    end
                end
				numBallObservations= numBallObservations-1;
				ballPos(numBallObservations,1)= blobs(num).blob.x(1);
				ballPos(numBallObservations,2)= blobs(num).blob.x(2);
				ballPos(numBallObservations,3)= blobs(num).blob.x(3);
				ballPos(numBallObservations,4)= (temp_time4-temp_time3)/1000/1000;
			end

		end% end data update this time.
            end

	end% end blob3 and blob1
	else
		%return 
    end
    end
        end
    end


	if numBallObservations < 1.0 * maxLen/2.0
        %return  %10
    end

    % judge flying direction
	ballFlyingDirection = 1; % from human towards robot

	for i = 1:numBallObservations-1
	
		if (ballPos(i,2)>= ballPos(i+1,2))  %Y direction
		
			break;
		end
	end
	if i ~= numBallObservations
	
		ballFlyingDirection = -1;  % from robot towards human

		for i = 1:numBallObservations-1
		
			if (ballPos(i,2) <= ballPos(i+1,2))
				return;
            end
		end
    end

    % find the lowest position
	zmin = 1000;
	zIndex = 1;
	upflag = TRUE;

	for i = 1:numBallObservations
	
	    if ballPos(i,3) < zmin		
			zmin = ballPos(i,3);
			zIndex = i;
		end
    end

    % select the predict mode
	if zIndex == 1 || zIndex == numBallObservations
	
		if ballFlyingDirection == 1 && ballPos(1,2) < -2.3
			PredictMode = 1;
        end
        
		if ballFlyingDirection == -1 && ballPos(1,2) > -2.1
			PredictMode = 1;
        end
        
		for i = 1:numBallObservations		
			if ballPos(i,3) > ballPos(i+1,3)
                upflag = false;
            end
        end
        
		if (ballFlyingDirection == 1 && ballPos(1,2) > -2.2 && upflag) % rebound
			PredictMode = 2;
        end
	end

	if zIndex ~= 1 && zIndex ~= numBallObservations
	
		if ballFlyingDirection= = 1 && ballPos(zIndex,2) <-2.2
            %return  % judge the landing position
        end
 		if ballFlyingDirection= = -1 && ballPos(zIndex,2)>-2.2) 
            %return 
        end

		%printf('\nget the rebound position\n');
		for temp1 = zIndex:numBallObservations
			for temp2 = 1:4
				ballPos(temp1+1-zIndex, temp2) = ballPos(temp1, temp2);
            end
        end
		numBallObservations = numBallObservations-zIndex+1;

		PredictMode = 2;
	end

% check the position number and determin the fitting length
	fitsize = 0; % fitting length (seen as the initial trajectory)
	if PredictMode == 1
	
		if numBallObservations < maxLen
			%return 
        end
		fitsize = numBallObservations;
    end
    
	if PredictMode == 2
	
		if numBallObservations < 1.0*maxLen/2.0
			%return 
        end
		fitsize = numBallObservations;
	end

% check flying direction again
	if (ballFlyingDirection == 1)
	
		for i = 1:numBallObservations-1
		
			if ballPos(i,2)>= ballPos(i+1,2)			
				%return 
			end
		end
	end
	if ballFlyingDirection == -1
	
		for i = 1:numBallObservations-1
		
			if ballPos(i,2) <= ballPos(i+1,2)			
				%return 
			end
		end
    end

    %2nd order polynomial: initial position/velocity
    px = zeros(4,1);
    py = zeros(4,1);
    pz = zeros(4,1);

	if polyfit(ballPos,fitsize,px,py,pz)
        % initial position
		back_No = fitsize/2.0;
		t_No = ballPos(back_No,4);
		ballInit.x(1) = px(1)*(t_No*t_No) + px(2)*t_No + px(3);
		ballInit.x(2) = py(1)*(t_No*t_No) + py(2)*t_No + py(3);
		ballInit.x(3) = pz(1)*(t_No*t_No) + pz(2)*t_No + pz(3);
		if abs(ballInit.x(1)-ballPos(back_No,1)) + abs(ballInit.x(2)-ballPos(back_No,2))+ ...
                abs(ballInit.x(3)-ballPos(back_No,3)) > 15.0/1000

			ballInit.x(1) = ballPos(back_No,1);
			ballInit.x(2) = ballPos(back_No,2);
			ballInit.x(3) = ballPos(back_No,3);
			%disp('change initial position');
        end

        % initial velocity
		ballInit.xd(1) = 2 * px(1) * t_No + px(2); % m/s
		ballInit.xd(2) = 2 * py(1) * t_No + py(2); % m/s
		ballInit.xd(3) = 2 * pz(1) * t_No + pz(2); % m/s

        % determine the past time
		ballTime = ballPos(numBallObservations,4) - ballPos(back_No,4); % predict from 'back_No'
    end

	if ballFlyingDirection * ballInit.xd(2)  < 0
		error('wrong direction!');
    end


