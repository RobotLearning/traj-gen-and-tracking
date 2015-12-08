classdef Phase < handle

      
   properties
      dt;
      z;
      zd;
      zdd;
      Tf;
   end
   
        
    
methods

    function obj = Phase(phase_dt, Tf)
        obj.dt = phase_dt;
        obj.Tf = Tf;
        
        obj = obj.derivatives;

    end

   function obj = derivatives(obj)
    % Phase and derivatives
 
        dt = obj.dt;
        
        z   = dt:dt:obj.Tf + 1*dt; 
        z = z';
        zD  = diff(z)./dt;
        zD  = [zD;  zD(end)];
        zDD = diff(zD)./dt;
        zDD = [zDD ; zDD(end)];

        obj.z   = z;
        obj.zd  = zD;
        obj.zdd = zDD;
   end
   
   
    
end %methods


end