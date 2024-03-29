function fode=TumorODE(chi,y,S,KOb)
fode=[y(2);
 S*y(1)/(KOb+y(1))-2/chi*y(2)];