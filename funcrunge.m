function Y = funcrunge(x,y,A,B,g_vet,u);

     y = [y(1);y(2);y(3);y(4)];
     u = [u(1);u(2)];
     Y = A*y - B*g_vet + B*u;
end
