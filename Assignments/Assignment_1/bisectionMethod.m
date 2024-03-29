function c = bisectionMethod(F,a,b,e)
c=(a+b)/2;
while abs(F(c))>e
    if F(c)<0&&F(a)<0
        a=c;
    else
        b=c;
    end
    c=(a+b)/2;
end