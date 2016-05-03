%% Finally disconnect the SL connection during clean up 

function disconnectFromSL(socket,address,context)

zmq.core.disconnect(socket, address);
zmq.core.close(socket);
zmq.core.ctx_shutdown(context);
zmq.core.ctx_term(context);

end