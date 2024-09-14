function Y = runge_kutta(x,y,h,A,B,g_vet,u);
%provavelmente u também
    k_1 = funcrunge(x,y,A,B,g_vet,u);
    k_2 = funcrunge(x+0.5*h,y+0.5*h*k_1,A,B,g_vet,u);
    k_3 = funcrunge((x+0.5*h),(y+0.5*h*k_2),A,B,g_vet,u);
    k_4 = funcrunge((x+h),(y+k_3*h),A,B,g_vet,u);
    Y = y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h; 
end

