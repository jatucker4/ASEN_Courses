function bin = Voltage2Bin(voltage, max_voltage, min_voltage, bits)
  nbins = 2^bits; % number of different levels

  % clip to max and min
  voltage = min([voltage, max_voltage]);
  voltage = max([voltage, min_voltage]);

  delta = (max_voltage - min_voltage)/(nbins); % resolution
  bin   = floor(voltage/delta); % make MATLAB do something 
                                % that resembles an unsigned 
                                % integer reading
end

