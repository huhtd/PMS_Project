% Generate a propagation movie

aviobj=avifile('mode_propagation.avi'); % The name of the output movie file
aviobj.fps=9; % Number of frames per second
aviobj.quality=100; % Quality of the movie
aviobj.compression='None'; % Compression algorithm (Unix compatible)

Px1=100;Px2=500;Py1=100;Py2=620;    % Size of figures
figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
set_figure_1;
Vc=Ez;
find_min_max;

fra=1;
for p=1:36 % loop to obtain 36 frames for a movie
p
    for i=1:Nn % loop over the grid to alter the field in time

        Vc(i)=Vc(i)*exp(imj*10*pi/180); % Phase shift (harmonic function in time)
   
    end %i

V=real(Vc);

% figure(2);
clf;
set_figure_1;
field_plot;
geometry_plot;
caxis([min_field max_field]);
colorbar('North');

M(fra)=getframe;
aviobj = addframe(aviobj,M(fra));
fra=fra+1;

end %p

movie(M,fra-1);

aviobj=close(aviobj);

return;
