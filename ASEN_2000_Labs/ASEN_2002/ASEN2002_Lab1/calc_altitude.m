function altitude = calc_altitude(density)
altitude_range = 35000:100:50000;
[T, a, P, rho] = atmoscoesa(altitude_range);
for i = 1:length(rho)
    altitude = altitude_range(i);
    if abs(rho(i) - density) < 0.0002
        break
    end
end
end