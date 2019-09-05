clear
clc
close all;
% C?lculo de las Series de Fourier

% --------------- VARIABLES GLOBALES ---------------
% iteraciones tabla
i=10;
% armonicos
arm=100;

% Constantes necesarias para el c?lculo
% periodo
T=12;
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



% -------- Calculo general de integrales --------
%
%
%
%

% Funcion por parte
tramo1=(0.5*t+2)*normalizador;
tramo2=(-0.2*t+0.3)*normalizador;
tramo3=(-0.3)*normalizador;
% Intervalos de tiempo de la funcion
t1 = -4;
t2 = -2;
t3 = 0;
t4 = 4;

% Calculo de las integrales basicas para valor medio
fx1 = int(tramo1,t1,t2);
fx2 = int(tramo2,t2,t3);
fx3 = int(tramo3,t3,t4);

% Calculo de las integrales para energ?a y potencia
fx1p = int(tramo1.^2,t1,t2);
fx2p = int(tramo2.^2,t2,t3);
fx3p = int(tramo3.^2,t3,t4);

energia = fx1p + fx2p + fx3p; % Se calcula la energia (suma de integrales)
potencia = (energia)*c; % c = 1/T. y Potencia = energia*1/T
valormedio = (fx1 + fx2 + fx3)*c; % c = 1/T

% print de resultados
fprintf("\n\nEnergia de la señal: %f | Potencia de la se?al: %f | Valor medio de la señal: %f \n\n", energia, potencia, valormedio);

% Calculo de las Integrales del coeficiente Fn de serie de Fourier Exponencial
Fx1 = int(tramo1*f_exponencial,t1,t2);
Fx2 = int(tramo2*f_exponencial,t2,t3);
Fx3 = int(tramo3*f_exponencial,t3,t4);

% Multiplicar por constante 1/T que es la constante global "c"
Fx1 = Fx1*c;
Fx2 = Fx2*c;
Fx3 = Fx3*c;

% Se suman para obtener el Fn
Fn = Fx1 + Fx2 + Fx3;

% -------- Calculo de tabla de coeficientes Fn + Coeficientes para graficos de espectro --------
%
%
%
%
fprintf("\n\n-------------------- TABLA DE COEFICIENTES FN --------------------\n\n");
coef_fn = zeros(1, 2*i+1); % creamos un vector de largo -i hasta i para los |fn^2|
coef_wn = zeros(1, 2*i+1); % creamos un vector de largo -i hasta i para los wn(rad/s)
coef_angle = zeros(1, 2*i+1); % creamos un vector de largo -i hasta i para los angle(rad)
for n=-i:-1
    [y, mf, wt, o] = coeficientesfn(Fn, n, w); % devuelve y=Fn, mf=|Fn^2|, wt=w*n, o=desfase
    coef_fn(n+i+1) = y; % El argumento n+i+1 sirve para llegar a la posicion 1 del array. Originalmente el N vale -i, si le sumamos i, vale 0, por lo tanto, falta sumarle un 1.
    coef_wn(n+i+1) = wt;
    coef_angle(n+i+1) = o;
    fprintf("Para n = %d | W(rad/s): %f | desf: %f | fn: %f%+fj | fn^2: %f\n", n, wt, o, real(y), imag(y), mf); % Se imprimen los valores
end

    % se evalua en un valor tendiendo a 0 para que no de error
    [y, mf] = coeficientesfn(Fn, 0.00001, w);
    coef_fn(i+1) = y;
    coef_wn(i+1) = 0;
    coef_angle(i+1) = 0;
    fprintf("Para n = 0 | W(rad/s): 0 | desf: 0 | fn: %f%+fj | fn^2: %f\n", real(y), imag(y), mf);
    
for n=1:i
    [y, mf, wt, o] = coeficientesfn(Fn, n, w);
    coef_fn(n+i+1) = y; % Como ya se escribió en los primeros i+1 posiciones del array, se lo sumamos a este argumento para comenzar una unidad despues.
    coef_wn(n+i+1) = wt;
    coef_angle(n+i+1) = o;
    fprintf("Para n = %d | W(rad/s): %f | desf: %f | fn: %f%+fj | fn^2: %f\n", n, wt, o, real(y), imag(y), mf);
end
%
%
% ------------------



% -------- Graficas de espectro de Fourier --------
%
%
%
%
% Grafico magnitud

figure; % Instancia donde se grafica, es decir, es la ventana que se abre al graficar.
hold on; % Mantener diferentes graficos en una figura.
title('Espectro magnitud de Fourier'); % Titulo del grafico
stem(coef_wn, abs(coef_fn.^2)); % grafico discreto stem del espectro de magnitud
ylabel('|Fn^2|'); % Titulo eje Y.
xlabel('w*n (rad/s)'); % Titulo eje X.
grid minor; % Cuadriculas en la grafica.

% Grafico fase

figure; % Instancia donde se grafica, es decir, es la ventana que se abre al graficar.
hold on; % Mantener diferentes graficos en una figura.
title('Espectro fase de Fourier'); % Titulo del grafico
stem(coef_wn, coef_angle); % grafico discreto stem del espectro de magnitud
ylabel('<º [-pi/2, pi/2] (rad)'); % Titulo eje Y.
xlabel('w*n (rad/s)'); % Titulo eje X.
grid minor; % Cuadriculas en la grafica.


% COMANDOS UTILES
% xlim([m n]) - Limita el eje X para valores de m a n
% ylim([m n]) - Limita el eje Y para valores de m a n
%
%
% ------------------



% -------- Calculo de grafica serie de Fourier --------
%
%
%
vector_t = -10:0.1:10;
fourier = zeros(1, length(vector_t));
fc = matlabFunction(Fn); % Transforma la ecuacion simbolica Fn calculada anteriormente a una función de matlab, sus argumentos equivalen a sus variables simbolicas.


for k=1:length(fourier) %Calculo de cada valor de la serie
    t = vector_t(k);
    suma = 0;
    for n=-arm:1:arm % Calculo de la sumatoria de la serie
        if n~=0
            suma = suma + fc(n)*exp(j*n*w*t); % Calculo del Fn y su evaluación con la exponencial, luego se añade al resultado de la sumatoria
        end
    end
    fourier(k) = suma + valormedio; % se le suma el valor medio que equivale a Fn con n=0.
end
figure;
hold on;
titulo = strcat("Serie de Fourier con " + arm + " armonicos");
title(titulo);
xlabel('t(s)');
ylabel('f(t)');
plot(vector_t, fourier);
grid minor;
%
%
% ------------------



% Declaracion e implementacion de funciones
% Funcion que devuelve distintos coeficientes
% y = Fn
% mf = |Fn^2|
% wt = w*n (rad/s)
% o = desfase [-pi/2, pi/2]
function [y, mf, wt, o] = coeficientesfn(Fn, n, w)
y=subs(Fn);
mf = abs(y^2);
wt = n*w;
o = atan(imag(y)/real(y));
end



