function  flecha(carga,xy)
lflecha=1.5;
sizecarga=size(carga,1);

xyarrow=xy;
for i=1:sizecarga  
    if carga(i,3)==1%% horizontal
        if carga(i,1)<0
            p0=[xyarrow(i,1)+lflecha,xyarrow(i,2)];
            p1=[xyarrow(i,1),xyarrow(i,2)];
            hArrow = drawArrow(p0,p1,'k');
        else
            p0=[xyarrow(i,1),xyarrow(i,2)];
            p1=[xyarrow(i,1)+lflecha,xyarrow(i,2)];
            hArrow = drawArrow(p0,p1,'k'); 
        end
       h = text(xyarrow(i,1)+lflecha,xyarrow(i,2), num2str(carga(i,1)));       
    end
    if carga(i,3)==2%% vertical
        if carga(i,1)<0
            p0=[xyarrow(i,1),xyarrow(i,2)];
            p1=[xyarrow(i,1),xyarrow(i,2)-lflecha];
            hArrow = drawArrow(p0,p1,'k');
        else
            p0=[xyarrow(i,1),xyarrow(i,2)-lflecha];
            p1=[xyarrow(i,1),xyarrow(i,2)];
            hArrow = drawArrow(p0,p1,'k'); 
        end
       h = text(xyarrow(i,1),xyarrow(i,2)-lflecha, num2str(carga(i,1)));       
    end
end
