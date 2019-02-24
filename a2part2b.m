
function I = a2part2b(PPS)
    %%%%%%%%%%%% Harmonic Wave Equation in 2D FD and Modes %%%%%%%%%%%%
    % 
    % Febuary 24th, 2019
    % Assignment 2
    % Patrobas Adewumi
    global C;
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                    % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458; % speed of light

    % Define area of region
    W = 10;    % width in y dir
    L = W*3/2; % length in x dir
    
    % Setting up the matrices for evaluation
    G = sparse(L*W,L*W);
    B = zeros(L*W,1);
    conduc = zeros(L,W);
    
    % Defining conductivity of the boxes (given area)
    s1 = 1;  % outiside boxes
    s2 = 0.01; % inside boxes

    % Define resistive region size
    res_L = L*1/4;
    res_W = W*2/5;

    % Centre point of given region
    mid_x = L/2;
    mid_y = W/2;
    
    % Investigating mesh density
    d_r = 1/PPS;  

    p_L = PPS*L;
    p_W = PPS*W;
    
    Smap = ones(L,W);
    for i = 1:1:p_L
        for j = 1:1:p_W
            n = j+(i-1)*p_W;
            nxm = j+(i-2)*p_W;
            nxp = j+i*p_W;
            nyp = j+1+ (i-1)*p_W;
            nym = j-1+ (i-1)*p_W;

            if(i == 1)
                conduc(i,j) = s1;
                G(n,:) = 0;
                G(n,n) = 1/(d_r^2);
                B(n) = 1;
                Smap(i,j) = s1;
            elseif(i == p_L)
                conduc(i,j) = s1;
                G(n,:) = 0;
                G(n,n) = 1/(d_r^2);
                B(n) = 1;
                Smap(i,j) = s1;
            elseif(j == 1)
                G(n,:) = 0;
                G(n,n) = -3/(d_r^2);
                if(i/PPS > mid_x - (res_L/2) && i/PPS < mid_x + (res_L/2))
                    conduc(i,j) = s1;
                    G(n,nxm) = s2/(d_r^2);
                    G(n,nxp) = s2/(d_r^2);
                    G(n,nyp) = s2/(d_r^2);
                    Smap(i,j) = s2;
                else
                    conduc(i,j) = s1;
                    G(n,nxm) = s1/(d_r^2);
                    G(n,nxp) = s1/(d_r^2);
                    G(n,nyp) = s1/(d_r^2);
                    Smap(i,j) = s1;
                end
            elseif(j == p_W)
                conduc(i,j) = s1;
                G(n,:) = 0;
                G(n,n) = -3/(d_r^2);
                if(i/PPS > mid_x - (res_L/2) && i/PPS < mid_x + (res_L/2))
                    G(n,nxm) = s2/(d_r^2);
                    G(n,nxp) = s2/(d_r^2);
                    G(n,nym) = s2/(d_r^2);
                    Smap(i,j) = s2;
                else
                    conduc(i,j) = s1;
                    G(n,nxm) = s1/(d_r^2);
                    G(n,nxp) = s1/(d_r^2);
                    G(n,nym) = s1/(d_r^2);
                    Smap(i,j) = s1;
                end
            else      
                conduc(i,j) = s1;
                G(n,:) = 0;
                G(n,n) = -4/(d_r^2);
                %  Setting my X and Y Boundaries                                | 
                if((i/PPS > mid_x -(res_L/2) && i/PPS < mid_x + (res_L/2)) && ...
                        (j/PPS > mid_y + (res_W/2) || j/PPS < mid_y - (res_W/2)))
           
                    conduc(i,j) = s1;
                    G(n,nxp) = s2/(d_r^2);
                    G(n,nxm) = s2/(d_r^2);
                    G(n,nyp) = s2/(d_r^2);
                    G(n,nym) = s2/(d_r^2);
                    Smap(i,j) = s2;
                else
                    conduc(i,j) = s1;
                    G(n,nxp) = s1/(d_r^2);
                    G(n,nxm) = s1/(d_r^2);
                    G(n,nyp) = s1/(d_r^2);
                    G(n,nym) = s1/(d_r^2);
                    Smap(i,j) = s1;
                end     
            end
        end
    end
%%
    V = G\B;

    Vmap = ones(p_L,p_W);
    for i = 1:1:p_L
        for j = 1:1:p_W
            n = j+(i-1)*p_W;
            Vmap(i,j) = V(n);
        end
    end

    E = gradient(Vmap);

    J = conduc.*Smap.*E;
    area = L*W;
    J_avg = sum(sum(J))/(p_L*p_W); % Current density
    I = J_avg/area;  % Current

end 