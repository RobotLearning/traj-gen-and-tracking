%% Decodes response from SL and outputs a structure
% containing arm positions and/or ball pos, camera id, time etc.
% TODO: get also camera status

function STR = decodeResponseFromSL(response)


offset = 1;
while offset < length(response)
    % get the command
    operationCode = response(offset);
    offset = offset + 1;
    switch operationCode
        case 0
            break; %means end of response
        case 3
            % observed arm positions
            vals = 4;
            numArmPos = typecast(response(offset:offset+vals-1),'uint32');
            offset = offset + vals;
            q = zeros(7,numArmPos);
            qd = zeros(7,numArmPos);
            t = zeros(1,numArmPos);
            for i = 1:numArmPos
                % get 15 doubles
                vals = 15*8;
                robotVals = typecast(response(offset:offset+vals-1),'double');
                q(:,i) = robotVals(1:7)';
                qd(:,i) = robotVals(8:14)';
                t(i) = robotVals(15);
                offset = offset + vals;
            end
            STR.robot.traj.q = q;
            STR.robot.traj.qd = qd;
            STR.robot.traj.time = t;
        case 4
            % get observed ball positions
            vals = 4;
            numBallPos = typecast(response(offset:offset+vals-1),'uint32');
            offset = offset + vals;
            ballPos = zeros(3,numBallPos);
            ballTime = zeros(1,numBallPos);
            camId = zeros(1,numBallPos);
            for i = 1:numBallPos
                % get 4 doubles
                vals = 4*8;
                ballVals = typecast(response(offset:offset+vals-1),'double');
                ballPos(:,i) = ballVals(1:3)';
                ballTime(i) = ballVals(4);
                offset = offset + vals;
                % get 2 integers
                vals = 2*4;
                ballVals = typecast(response(offset:offset+vals-1),'uint32');
                camId(i) = ballVals(1);
                % 2nd integer is for padding
                offset = offset + vals;
            end   
            STR.ball.pos = ballPos;
            STR.ball.time = ballTime;
            STR.ball.cam.id = camId;
        case 8 
            % get robot time
            vals = 8;
            tRobot = typecast(response(offset:offset+vals-1),'double');
            offset = offset + vals;
            STR.robot.time = tRobot;
        otherwise
            error('Operation code not recognized!');
    end
            
end