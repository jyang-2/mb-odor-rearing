pentanol_conc = [0, -3, -4, -5, -6];
pentanol_pins = containers.Map(pentanol_conc,[nan, nan, nan, nan, nan]);
pentanol_pins(0) = 42;  % pfo
pentanol_pins(-6) = 41;
pentanol_pins(-5) = 40;
pentanol_pins(-4) = 39;
pentanol_pins(-3) = 38;

hexanol_conc = [0, -3, -4, -5, -6];
hexanol_pins = containers.Map(hexanol_conc,[nan, nan, nan, nan, nan]);
hexanol_pins(0) = 49;    % pfo
hexanol_pins(-6) = 48;
hexanol_pins(-5) = 47;
hexanol_pins(-4) = 46;
hexanol_pins(-3) = 45;
%%
conc1_pins = containers.Map(pentanol_pins.values(), pentanol_pins.keys());
conc2_pins = containers.Map(hexanol_pins.values(), hexanol_pins.keys());

conc = [0,0,0,-6,-6,-6,-5,-5,-5,-4,-4,-4,-3,-3,-3,-6,-6,-6,-5,-5,-5,-4,-4,-4,-3,-3,-3];

channelA = [42,42,42,41,42,41,40,42,40,42,39,39,38,38,42,42,41,41,40,42,40,42,39,39,38,38,42];
channelB = [49,49,49,48,48,49,49,47,47,46,49,46,45,49,45,48,49,48,47,47,49,46,49,46,49,45,45];


conc1 = arrayfun(@(x) conc1_pins(x), channelA);
conc2 = arrayfun(@(x) conc2_pins(x), channelB);
