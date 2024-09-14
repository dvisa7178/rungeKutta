function [theta1, theta2] = calcAngulos(x, y, a1, a2)
    cos_theta2 = (x^2 + y^2 - a1^2 - a2^2) / (2 * a1 * a2);
    theta2 = acos(cos_theta2);
    theta1 = atan2(y, x) - atan2(a2 * sin(theta2), a1 + a2 * cos(theta2));
end
