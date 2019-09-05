clear
clc
% C?lculo de las Series de Fourier

% --------------- VARIABLES GLOBALES ---------------
% iteraciones tabla
i=10;
% armonicos
arm=100;

% Constantes necesarias para el c?lculo
% periodo
T=16;
% frecuencia
w=2*pi/T;

% --------------- NO TOCAR SI NO SABE ---------------
% Variables simbolicas
syms t;
syms n;
% Constantes globales
c = 1/T;
normalizador = t/t; % utilizado para que matlab no arroje error en algunas integrales.
f_exponencial = exp(-j*n*t*w);


% --------------- INICIO DE CODIGO LOGICO ---------------
% Funci?n por parte
tramo1=(11)*normalizador;
tramo2=(2.7183)*normalizador;
tramo3=(1)*normalizador;
% Intervalos de tiempo
t1 = -8;
t2 = -7;
t3 = -3;
t4 = 5;

% C?lculo de las Integrales b?sicas
% para valor medio
fx1 = int(tramo1,t1,t2);
fx2 = int(tramo2,t2,t3);
fx3 = int(tramo3,t3,t4);
% para energ?a y potencia
fx1p = int(tramo1.^2,t1,t2);
fx2p = int(tramo2.^2,t2,t3);
fx3p = int(tramo3.^2,t3,t4);

energia = fx1p + fx2p + fx3p;
potencia = (energia)*c;
valormedio = (fx1 + fx2 + fx3)*c;
fprintf("\n\nEnergia de la se?al: %f | Potencia de la se?al: %f | Valor medio de la se?al: %f \n\n", energia, potencia, valormedio);

% C?lculo de las Integrales de la serie de Fourier
Fx1 = int(tramo1*f_exponencial,t1,t2);
Fx2 = int(tramo2*f_exponencial,t2,t3);
Fx3 = int(tramo3*f_exponencial,t3,t4);

% Multiplicar por constante 1/T que es la constante global "c"
Fx1 = Fx1*c;
Fx2 = Fx2*c;
Fx3 = Fx3*c;

% Fn final

Fn = Fx1 + Fx2 + Fx3;

% C?lculo de tabla de coeficientes Fn
fprintf("\n\n-------------------- TABLA DE COEFICIENTES FN --------------------\n\n");
for n=-i:-1
    [y, mf, wt, o] = coeficientesfn(Fn, n, w);
    fprintf("Para n = %d | W(rad/s): %f | desf: %f | fn: %f%+fj | fn^2: %f\n", n, wt, o, real(y), imag(y), mf);
end
    [y, mf] = coeficientesfn(Fn, 0.00001, w);
    fprintf("Para n = 0 | W(rad/s): 0 | desf: 0 | fn: %f%+fj | fn^2: %f\n", real(y), imag(y), mf);
for n=1:i
    [y, mf, wt, o] = coeficientesfn(Fn, n, w);
    fprintf("Para n = %d | W(rad/s): %f | desf: %f | fn: %f%+fj | fn^2: %f\n", n, wt, o, real(y), imag(y), mf);
end


% C?lculo de grafica serie de Fourier
vector_t = -27:0.1:27;
fourier = zeros(1, length(vector_t));
fc = matlabFunction(Fn);
for k=1:length(fourier)
    t = vector_t(k);
    suma = 0;
    for n=-arm:1:arm
        if n~=0
            suma = suma + fc(n)*exp(j*n*w*t);
        end
    end
    fourier(k) = suma + valormedio;
end
plot(vector_t, fourier);



% Declaracion e implementacion de funciones
function [y, mf, wt, o] = coeficientesfn(Fn, n, w)
y=subs(Fn);
mf = abs(y^2);
wt = n*w;
o = atan(imag(y)/real(y));
end



