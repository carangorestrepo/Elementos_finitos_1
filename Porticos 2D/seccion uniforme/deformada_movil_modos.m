datos=100;
def1=zeros(elementos(1,1)*6,1);
recorrido=1:elementos(1,1):elementos(1,1)*6;
def=zeros(nudos*6,1);
xyed=zeros(nudos*4,4);
modo=2;
escala=1;
w=Omega(modo)*escala;

def(d,1)=Phi(:,modo)*10;
esch=2;
escv=10;
%h = zeros(elemetos,1);
t=0:0.005:2;
xmax=max(abs(def))+max(max(abs([xye(:,1),xye(:,3)])));
xmin=max(abs(def))+min(min(abs([xye(:,1),xye(:,3)])));
ymax=max(abs(def))+max(max(abs([xye(:,2),xye(:,4)])));
ymin=max(abs(def))+min(min(abs([xye(:,2),xye(:,4)])));

ymin=0;
grid on;
myVideo = VideoWriter('myVideoFile-4'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)

for r=1:datos
    hold on;
    for e=1:elementos
        Xi=def(GLe(e,1))*cos(w*t(r));
        Yi=def(GLe(e,2))*cos(w*t(r));
        Xf=def(GLe(e,4))*cos(w*t(r));
        Yf=def(GLe(e,5))*cos(w*t(r));
        xyed(e,:) = [xye(e,1)+Xi,xye(e,2)+Yi,xye(e,3)+Xf,xye(e,4)+Yf];
        grid on
        h(e)=plot(xyed(e,[1,3]),xyed(e,[2,4]),'--k','LineWidth', 2);hold on;
        axis equal
    end
    xlim([-xmin,xmax]);
    ylim([ymin,ymax]);

    lx=num2str(T(modo));
    by=num2str(w/(2*pi()));
    nombre=double('Modo de vibraciòn');
    nombre1=double('/');
    nombre=[nombre,32,num2str(modo),32,double('f='),by,double('(Hz)'),nombre1,double('T='),lx,double('(s)')];
    xlabel('x(m)')
    ylabel('y(m)')
    zlabel('z(m)')
    
    title(char(nombre));
    pause(0.1) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    
    drawnow;
    %pause(0.005);
    delete(h)
    hold off;
    %cla(h)

end
close(myVideo)