%% Test communication to SL via matlab-zmq

host = 'localhost'; 
port = '7646';
address = sprintf('tcp://%s:%s',host,port);
context = zmq.core.ctx_new();
socket  = zmq.core.socket(context, 'ZMQ_REQ');
zmq.core.connect(socket, address);

msg = [uint8(3), uint32(1), uint8(0)];
data = typecast(msg, 'uint8'); 
zmq.core.send(socket, data);
response = zmq.core.recv(socket);

STR = decodeResponseFromSL(response);

zmq.core.disconnect(socket, address);
zmq.core.close(socket);
zmq.core.ctx_shutdown(context);
zmq.core.ctx_term(context);