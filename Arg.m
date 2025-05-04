function [arg] = Arg(sinT,cosT)
% Calculate the argument related to the time-delay τ
len=length(sinT);
arg=linspace(0,0,len);
n=0;
for i=1:len
    if i>=2
    if atan2(sinT(i-1),cosT(i-1))<-pi/2 && atan2(sinT(i),cosT(i))>pi/2
        n=n-1;
    elseif atan2(sinT(i-1),cosT(i-1))>pi/2 && atan2(sinT(i),cosT(i))<-pi/2
        n=n+1;
    end
    end
    arg(i)=atan2(sinT(i),cosT(i))+n*2*pi;
end
end