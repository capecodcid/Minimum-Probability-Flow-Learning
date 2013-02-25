function kind = threepointkind(cix, nN, nT)
    ntt = size(cix,1);
    kind = zeros(ntt,1);
    
    ixx = 1;
    index = zeros(nT * (nN^2-1),3);
    for t=1:nT,
        for x=1:nN,
           for y=1:nN,
               if (x == (nN+1)/2 & y==(nN+1)/2) continue; end;
               index(ixx,1)=x;
               index(ixx,2)=y;
               index(ixx,3)=t;
               ixx=ixx+1;
           end
        end
    end
           
    
    for i=1:ntt,
        ids = cix(i,:);
        c1 = index(ids(1),:);
        c2 = index(ids(2),:);
        c3 = index(ids(3),:);
        x1=c1(1);y1=c1(2);t1=c1(3);
        x2=c2(1);y2=c2(2);t2=c2(3);
        x3=c3(1);y3=c3(2);t3=c3(3);
        
        %abs(x1-x2) + abs(x1-x3) + abs(y1-y2) + abs(y1-y3) + abs(x2-x3) + abs(y2-y3)
        
        if (t1==t2 & t2==t3) kind(i) = 1;
        elseif ( abs(x1-x2) + abs(x1-x3) + abs(y1-y2) + abs(y1-y3) + abs(x2-x3) + abs(y2-y3) == 2)
            kind(i) = 2;
        else kind(i) = 0;
        end
    end
end
