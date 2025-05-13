

modo=1;
w=Omega(modo);
datos=100;
t=linspace(0,2,datos);
def=zeros(nudos*6,1);
xyed=zeros(nudos*6,6);
escx=2;
escy=2;
escz=10;

def(GLKM,1)=Phi(:,modo)*10;

xmax=max(abs(def))+max(abs(xye(:,4)));
xmin=max(abs(def))+min(abs(xye(:,1)));
ymax=max(abs(def))+max(abs(xye(:,2)));
ymin=max(abs(def))+min(abs(xye(:,5)));

myVideo = VideoWriter('myVideoFile-4'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)

for r=1:datos
    for e=1:elementos(1,1)
        Xi=def(GLe(e,1))*cos(w*t(r));
        Zi=def(GLe(e,2))*cos(w*t(r));
        Yi=def(GLe(e,3))*cos(w*t(r));
        Xf=def(GLe(e,7))*cos(w*t(r));
        Zf=def(GLe(e,8))*cos(w*t(r));
        Yf=def(GLe(e,9))*cos(w*t(r));
        xyed(e,:) = [xye(e,1)+Xi,xye(e,2)+Yi,xye(e,3)+Zi,xye(e,4)+Xf,xye(e,5)+Yf,xye(e,6)+Zf];
        grid on
        h(i)=plot3(xyed(e,[1,4]),xyed(e,[2,5]),xyed(e,[3,6]),'--k','LineWidth', 2);hold on;
        axis equal
    end
    xlim([-xmin,xmax]);
    ylim([-ymin,ymax]);
    
    

    
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    
    drawnow;
    %pause(0.005);
    hold off;
    %cla(h)
    delete(h);
end


close(myVideo)