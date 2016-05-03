%% Observe ball from SL

msg = [uint8(4), typecast(uint32(1),'uint8'),uint8(0)];
data = typecast(msg,'uint8');
zmq.core.send(socket,data);
response = zmq.core.recv(socket,bufferLength);
STR = decodeResponseFromSL(response);
ballObs = STR.ball.pos;
ballTime = STR.ball.time;
ballCam = STR.ball.cam.id;