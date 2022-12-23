function Ez=te_even_mode(y,E0,Order,d)

% y     - Variable, Ez=Ez(y)
% Order - Order of the mode
% d     - Thicknes of the waveguide
% E0    - Amplitude of the signal

Ez=E0*cos((2*Order+1)*pi/d*y);

end