function Y = funcrunge(x,y,A,B,g);

     y = [y(1);y(2);y(3);y(4)];
     u = [0.0;0.0];
     Y = A*y - B*g + B*u;
end
