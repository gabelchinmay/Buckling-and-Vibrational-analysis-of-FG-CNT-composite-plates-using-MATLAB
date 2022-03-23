clc
clear all
close all
syms z
syms x y a b r s

% ---------------------------------------- SETTING THE PARAMETERS AND [ABD]-------------------------------------------
%defining the dimensions of the plate (in m)
len = 1; %length
br = len; % breadth
%th = 2e-3; % thickness
th = len/50;
% defining the coordinates of individual elements
x1 = 0;
x2 = 2*a;
x3 = 2*a;
x4 = 0;
y1 = 0;
y2 = 0;
y3 = 2*b;
y4 = 2*b;

%defining the properties of steel (ASTM A36)
% E = 200e9; % Youngs Modulus in Pa
% G = 79.3e9; % Shear Modulus in Pa
% mu = 0.26; % possion's ratio
% rho =  7850; %density in kg/m^3

%Defining the properties of the matrix ( of the composite)PmPv
mu_m = 0.34;
rho_m = 1150;
E_m = 2.1e9;
G_m = E_m/(2*(1+mu_m));

%Defining the properties of the fibers (CNT)
mu_cnt = 0.175;
rho_cnt = 1400;
G12_cnt = 1.9445e12;
E11_cnt = 5.6466e12;
E22_cnt = 7.08e12;

% CNT effeciency parameter
V_cnt_prime = 0.17;
n1 = 0.149;
n2 = 1.381;
n3 = 1.381;

%Functionally graded distribution
V_cnt = V_cnt_prime*2*((2*abs(z)/th));
E11 = n1*V_cnt*E11_cnt + (1 - V_cnt)*E_m;
E22 = n2/((V_cnt/E22_cnt) + ((1 - V_cnt)/E_m));
G12 = n3/((V_cnt/G12_cnt) + ((1 - V_cnt)/G_m));
mu12 = V_cnt_prime*mu_cnt +  (1 - V_cnt)*mu_m;
rho = V_cnt*rho_cnt +(1 - V_cnt)*rho_m;



% defining the property matrix
%A = ([1/E -mu/E 0; -mu/E 1/E 0; 0 0 (1+mu)/E]^(-1));
S=[1/E11 -mu12/E11 0; -mu12/E11 1/E22 0; 0 0 1/G12];
Q=inv(S);
A=Q;
B = -z.*A;
D = (z^2).*A;
A = int(A,z,-th/2,th/2);
B = int(B,z,-th/2,th/2);
D = int(D,z,-th/2,th/2);
ABD = [A B; B D]; % the ABD matrix which contains property terms and 'z' terms
%ABD = vpa(simplify(ABD));

% ---------------------------------------- SHAPE FUNCTIONS -------------------------------------------
%defining the shape functions ( for u and v)
F = [1 x y x*y]*([1 0 0 0;1 2*a 0 0;1 2*a 2*b 4*a*b;1 0 2*b 0]^-1);
F = subs(F,x, a*r + a);
F = subs(F,y,b*s+b);
F = simplify(F); % this contains the shape functions for both u and v

% Denifing the shape functions (for w)
G = [1 x y x^2 x*y y^2 x^3 (x^2)*y x*(y^2) y^3 (x^3)*y x*(y^3)]; %this is for w
G1 = diff(G,y); % This is for the theta x
G2 = -diff(G,x); % This is for theta y
G_final = [G;G1;G2]; 
H = zeros(0,0);
for i = 1:4
    G_f = G_final;
    if i == 1
        G_f = subs(G_f,x,x1);
        G_f = subs(G_f,y,y1);
        H = [H;G_f];
    elseif i == 2
        G_f = subs(G_f,x,x2);
        G_f = subs(G_f,y,y2);
        H = [H;G_f];
    elseif i == 3
        G_f = subs(G_f,x,x3);
        G_f = subs(G_f,y,y3);
        H = [H;G_f];
    elseif i == 4
       G_f = subs(G_f,x,x4);
       G_f = subs(G_f,y,y4);
        H = [H;G_f];
    end
end
I = G*(H^-1);
I = subs(I,x, a*r + a);
I = subs(I,y,b*s+b);
I = simplify(I); % contains shape functions for w

% ---------------------------------------- FEA-------------------------------------------

% forming the matrix having Ns
J = sym(zeros(6,20)); %  initializing matrix having various derative of the shape functions 

for i = 1:6
    if i == 1
        k = 1;
        for j = 1:5:20
            J(i,j) = diff(F(k),r)/a;
            k = k + 1;
        end
    elseif i == 2
        k = 1;
        for j = 2:5:20
            J(i,j) = diff(F(k),s)/b;
            k = k+1;
        end
    elseif i == 3
        k = 1;
        for j = 1:5:20
            J(i,j) = diff(F(k),s)/b;
            J(i,j+1) = diff(F(k),r)/a;
            k = k + 1;
        end
    elseif i == 4
        k = 1;
        for j = 3:5:20
            J(i,j) = diff(diff(I(k),r),r)/(a^2);
            J(i,j+1) = diff(diff(I(k+1),r),r)/(a^2);
            J(i,j+2) = diff(diff(I(k+2),r),r)/(a^2);
            k = k + 3;
        end
    elseif i == 5
        k = 1;
        for j = 3:5:20
            J(i,j) = diff(diff(I(k),s),s)/(b^2);
            J(i,j+1) = diff(diff(I(k+1),s),s)/(b^2);
            J(i,j+2) = diff(diff(I(k+2),s),s)/(b^2);
            k = k + 3;
        end
    elseif i == 6
        k = 1;
        for j = 3:5:20
            J(i,j) = 2*diff(diff(I(k),r),s)/(a*b);
            J(i,j+1) = 2*diff(diff(I(k+1),r),s)/(a*b);
            J(i,j+2) = 2*diff(diff(I(k+2),r),s)/(a*b);
            k = k + 3;
        end
    end
end
J = simplify(J); % matrix having various derative of the shape functions 
n = input('Enter the value of n:'); % no. of nodes along x and y directions

% ---------------------------------------- STIFFNESS MATRIX-------------------------------------------

%getting the local stiffness matrix
stiffness_matrix = transpose(J)*ABD*J;
stiffness_matrix = int(int(stiffness_matrix,r,-1,1),s,-1,1); % integrating
stiffness_matrix = subs(stiffness_matrix,a,len/(2*(n-1))); % substituting the value of a
stiffness_matrix = subs(stiffness_matrix,b,br/(2*(n-1))); % substituting the value of b
stiffness_matrix_local = double(vpa(stiffness_matrix)); % simplifying

% discreatization of the nodes
Node_mat = zeros(n,n); % Initializing the matrix containing the nodes
Coll_mat = zeros(0,0); % Initializing the matrix containing the node info about all the elements
count  = 1;
for i = 1:n
    for j = n:-1:1
        Node_mat(j,i) = count; % matrix containing the nodes
        count = count + 1;
    end
end

% getting the node info of every element
for i = 1:n-1
    for j = n:-1:2
        c_mat = [Node_mat(j,i) Node_mat(j,i+1) Node_mat(j-1,i+1) Node_mat(j-1,i)];
        Coll_mat = [Coll_mat; c_mat]; %matrix containing the node info about all the elements
    end
end

% making the global stiffness matrix
stiffness_matrix_global = zeros(5*n*n,5*n*n); % initializing the global stiffness matrix

for i = 1:(n-1)^2
    m = Coll_mat(i,:);
    for j = 1:4
        ele1 =  m(j);
        for k = 1:4
            ele2 = m(k);
            stiffness_matrix_global(5*ele1-4:5*ele1,5*ele2-4:5*ele2) = stiffness_matrix_local(j*5-4:j*5,k*5-4:5*k) + stiffness_matrix_global(5*ele1-4:5*ele1,5*ele2-4:5*ele2);
        end
    end
end


% ---------------------------------------- MASS MATRIX-------------------------------------------

mass_1 = rho*[1 0 0 -z 0;0 1 0 0 -z;0 0 1 0 0;-z 0 0 z^2 0;0 -z 0 0 z^2];


% forming the matrix having Ns
J_mass = sym(zeros(5,20)); %  initializing matrix having various derative of the shape functions 

for i = 1:5
    if i == 1
        k = 1;
        for j = 1:5:20
            J_mass(i,j) = F(k);
            k = k + 1;
        end
    elseif i == 2
        k = 1;
        for j = 2:5:20
            J_mass(i,j) = F(k);
            k = k+1;
        end
    
    elseif i == 3
        k = 1;
        for j = 3:5:20
            J_mass(i,j) = I(k);
            J_mass(i,j+1) = I(k+1);
            J_mass(i,j+2) = I(k+2);
            k = k + 3;
        end
    elseif i == 4
        k = 1;
        for j = 3:5:20
            J_mass(i,j) = diff(I(k),r)/a;
            J_mass(i,j+1) = diff(I(k+1),r)/a;
            J_mass(i,j+2) = diff(I(k+2),r)/a;
            k = k + 3;
        end
    elseif i == 5
        k = 1;
        for j = 3:5:20
            J_mass(i,j) = diff(I(k),s)/(b);
            J_mass(i,j+1) = diff(I(k+1),s)/(b);
            J_mass(i,j+2) = diff(I(k+2),s)/(b);
            k = k + 3;
        end
    end
end
J_mass = simplify(J_mass); % matrix having various derative of the shape functions 
mass_1 = simplify(int(mass_1,z,-th/2,th/2)); % integrating the mass_1 matrix wrt z


%getting the local mass matrix
mass_matrix = transpose(J_mass)*mass_1*J_mass;
mass_matrix = int(int(mass_matrix,r,-1,1),s,-1,1); % integrating
mass_matrix = subs(mass_matrix,a,len/(2*(n-1))); % substituting the value of a
mass_matrix = subs(mass_matrix,b,br/(2*(n-1))); % substituting the value of b
mass_matrix_local = double(vpa(mass_matrix)); % simplifying


% making the global mass matrix
mass_matrix_global = zeros(5*n*n,5*n*n); % initializing the global mass matrix

for i = 1:(n-1)^2
    m = Coll_mat(i,:);
    for j = 1:4
        ele1 =  m(j);
        for k = 1:4
            ele2 = m(k);
            mass_matrix_global(5*ele1-4:5*ele1,5*ele2-4:5*ele2) = mass_matrix_local(j*5-4:j*5,k*5-4:5*k) + mass_matrix_global(5*ele1-4:5*ele1,5*ele2-4:5*ele2);
        end
    end
end

% ---------------------------IMPOSING BOUNDARY CONDITIONS-------------------
%  stiffness_matrix_global(1:n*5,:) = [];
%  stiffness_matrix_global(:,1:n*5) = [];
%  mass_matrix_global(1:n*5,:) = [];
%  mass_matrix_global(:,1:n*5) = [];
stiffness_matrix_global = stiffness_matrix_global(n*5+1:(n*n*5 -n*5),n*5+1:(n*n*5 -n*5));
mass_matrix_global = mass_matrix_global(n*5+1:(n*n*5 -n*5),n*5+1:(n*n*5 -n*5));

for i = 1 : (n-2)
        stiffness_matrix_global(:, 5*(n-1)*i+1 : 5*(n-1)*i+5) = [];
        mass_matrix_global(: ,5*(n-1)*i+1 : 5*(n-1)*i+5) = [];
end

for i = 1 : (n-2)
        stiffness_matrix_global(5*(n-1)*i+1 : 5*(n-1)*i+5,:) = [];
        mass_matrix_global(5*(n-1)*i+1 : 5*(n-1)*i+5,:) = [];
end

for i = 0 : (n-2)-1
        stiffness_matrix_global(:,5*(n-2)*i+1 : 5*(n-2)*i+5) = [];
        mass_matrix_global(:,5*(n-2)*i+1 : 5*(n-2)*i+5) = [];
end

for i = 0 : (n-2)-1
    stiffness_matrix_global(5*(n-2)*i+1 : 5*(n-2)*i+5,:) = [];
    mass_matrix_global(5*(n-2)*i+1 : 5*(n-2)*i+5,:) = [];
end

[eigvec,eigval]=eig(stiffness_matrix_global,mass_matrix_global);
lambda = (sqrt(diag(eigval)));%/(2*pi);
lambda_nd = lambda*((len^2/th)*sqrt(rho_m/E_m));

% siz = size(stiffness_matrix_global);
% ms = eigvec(3:5:siz(1),1);
% a1 = max(abs(ms));
% ms = (1/a1)*ms;
% c2 =0;
% for k =1:11
%     mss(1:n,k) = ms((c2*n)+1:(c2+1)*n);
%     c2=c2+1;
% end
% modeshape=[zeros(n,1) mss];
% surf(modeshape)