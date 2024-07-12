function Y = runge_kutta(x,y,h);
    k_1 = funcrunge(x,y);
    k_2 = funcrunge(x+0.5*h,y+0.5*h*k_1);
    k_3 = funcrunge((x+0.5*h),(y+0.5*h*k_2));
    k_4 = funcrunge((x+h),(y+k_3*h));
    Y = y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h; 
    
  
