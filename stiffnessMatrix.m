function ke = stiffnessMatrix(a,b,t,E,nu)

    % % % Input:
    % a = Element length in x-direction
    % b = Element length in y-direction
    % t = Thickness
    % E = Young's modulus
    % nu = Poisson's ratio
    % ey = Number of elements in y-direction
    % % % Output:
    % ke = Element stiffness matrix (8x8)
    % % %
    
t1 = -nu ^ 2 + 1;
t2 = 1 - nu;
t3 = 1 / b;
t1 = 1 / t1;
t4 = 1 / a;
t5 = t4 ^ 2;
t6 = t3 ^ 2;
t7 = t5 / 0.16e2;
t8 = t2 / 0.32e2;
t9 = t8 * t6;
t1 = E * t1;
t10 = t1 * (t9 + t7);
t11 = (a * t2 * t3) / 0.24e2;
t12 = b * t4;
t13 = t12 / 0.12e2;
t14 = t1 * (t13 + t11);
t15 = t * (0.4e1 * t10 * a * b + t14);
t16 = nu / 0.16e2;
t4 = t1 * t4 * t3;
t17 = t4 * (t8 + t16);
t7 = t1 * (t9 - t7);
t9 = t1 * (-t13 - t11);
t11 = t * (0.4e1 * t7 * a * b + t9);
t4 = (t4 * (-t8 + t16));
t10 = t * (-0.4e1 * t10 * a * b + t14);
t7 = t * (-0.4e1 * t7 * a * b + t9);
t9 = 0.4e1 * t;
t13 = t9 * t17 * a * b;
t14 = t9 * t4 * a * b;
t16 = -t9 * t17 * a * b;
t4 = -(t9 * t4 * a * b);
t6 = t6 / 0.16e2;
t5 = t8 * t5;
t8 = t1 * (t5 + t6);
t2 = (t12 * t2 / 0.24e2);
t3 = ((a * t3) / 0.12e2);
t9 = (t1 * (t3 + t2));
t12 = t * (0.4e1 * t8 * a * b + t9);
t5 = t1 * (-t5 + t6);
t1 = t1 * (-t3 - t2);
t2 = (t * (0.4e1 * t5 * a * b + t1));
t3 = (t * (-0.4e1 * t8 * a * b + t9));
t1 = (t * (-0.4e1 * t5 * a * b + t1));
ke = [t15 t13 t11 t14 t10 t16 t7 t4; t13 t12 t4 t2 t16 t3 t14 t1; t11 t4 t15 t16 t7 t14 t10 t13; t14 t2 t16 t12 t4 t1 t13 t3; t10 t16 t7 t4 t15 t13 t11 t14; t16 t3 t14 t1 t13 t12 t4 t2; t7 t14 t10 t13 t11 t4 t15 t16; t4 t1 t13 t3 t14 t2 t16 t12;];

end