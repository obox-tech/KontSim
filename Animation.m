function Animation(y,t)

%Animation
figure('position',[0 0 1000 700]);
axis(gca,'equal');
axis([-0.0 70.0 0.0 5.0]);
ax = gca;
ax.Units = 'pixels';
pos = ax.Position;
marg = 30;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
grid on;


for i=1:length(t)
    %Position Kart
    Kart_pos=[y(i,1) 0];
    %Position Gelenk
    Gelenk_pos=Kart_pos+[0.5*sin(y(i,2)) 0.5*cos(y(i,2))];
    %Position Ende
    Ende_pos=Gelenk_pos+[0.7*sin(y(i,3)) 0.7*cos(y(i,3))];
    
    %Kart als Rechteck
    K_point=rectangle('Position',[Kart_pos-[0.1 0.025] 0.2 0.05]);
    %Gelenk als Kreis
    G_point=viscircles(Gelenk_pos,0.01);
    %Stab 1
    stab1=line([Kart_pos(1) Gelenk_pos(1)],[Kart_pos(2) Gelenk_pos(2)]);
    %Stab 2
    stab2=line([Gelenk_pos(1) Ende_pos(1)],[Gelenk_pos(2) Ende_pos(2)]);
    %Ende als Kreis
    E_point=viscircles(Ende_pos,0.01);
    
    %Time interval to update the plot
    %pause(0.001);
    
    mov(i)=getframe(gcf);
    
    %Delete previous objects
    if i<length(t)
        delete(K_point);
        delete(G_point);
        delete(E_point);
        delete(stab1);
        delete(stab2);
    end
          
end

%movie(mov,1,200);
%movie2avi(mov,'testanimation3.avi')
v=VideoWriter('testanimation4');
v.FrameRate = 50;
v.Quality = 100;
open(v);
writeVideo(v,mov);
close(v);


end