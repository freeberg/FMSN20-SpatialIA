sigma = 0.35;
epsilon = 0.25;

p_x8 = (1-sigma)^7 * sigma^1 * (1-epsilon)^4 * epsilon^3
p_x0 = (1-sigma)^5 * sigma^2 * (1-epsilon)^3 * epsilon^5
p_xB = (1-sigma)^5 * sigma^0 * (1-epsilon)^5 * epsilon^5

p_x = p_xB * 0.25 + p_x0 * 0.4 + p_x8 * 0.35

p_Bx = p_xB * 0.25 / p_x
p_0x = p_x0 * 0.40 / p_x
p_8x = p_x8 * 0.35 / p_x