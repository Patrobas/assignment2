function I = a2part2d(sigma)
        
    %%%%%%%%%%%% Harmonic Wave Equation in 2D FD and Modes %%%%%%%%%%%%
    % 
    % Febuary 24th, 2019
    % Assignment 2
    % Patrobas Adewumi

        global C
        C.q_0 = 1.60217653e-19;             % electron charge
        C.hb = 1.054571596e-34;             % Dirac constant
        C.h = C.hb * 2 * pi;                % Planck constant
        C.m_0 = 9.10938215e-31;             % electron mass
        C.kb = 1.3806504e-23;               % Boltzmann constant
        C.eps_0 = 8.854187817e-12;          % vacuum permittivity
        C.mu_0 = 1.2566370614e-6;           % vacuum permeability
        C.c = 299792458;                    % speed of light
        C.g = 9.80665;


         % Define area of region
    W = 50;    % width in y dir
    L = W*3/2; % length in x dir
 
      % Centre point of given region
    mid_x = L/2;
    mid_y = W/2;
    
  % Setting up the matrices for evaluation
    G = zeros(L*W,L*W);
    B = zeros(L*W,1);

   % Defining conductivity of the boxes (given area)
    s1 = 1;
    s2 = sigma;

    % Define resistive region size
    res_L = L*1/4;
    res_W = W*2/5;

    Smap = zeros(L,W);
    for i = 1:1:L
        for j = 1:1:W
            n = j+(i-1)*W;
            nxm = j+(i-2)*W;
            nxp = j+i*W;
            nyp = j+1+ (i-1)*W;
            nym = j-1+ (i-1)*W;

            if(i == 1)
                G(n,:) = 0;
                G(n,n) = 1;
                B(n) = 1;
                Smap(i,j) = s1;
            elseif(i == L)
                G(n,:) = 0;
                G(n,n) = 1;
                B(n) = 0;
                Smap(i,j) = s1;
            elseif(j == 1)
                G(n,:) = 0;
                G(n,n) = -3;
                if(i > mid_x - (res_L/2) && i < mid_x + (res_L/2))
                    G(n,nxm) = s2;
                    G(n,nxp) = s2;
                    G(n,nyp) = s2;
                    Smap(i,j) = s2;
                else
                    G(n,nxm) = s1;
                    G(n,nxp) = s1;
                    G(n,nyp) = s1;
                    Smap(i,j) = s1;
                end
            elseif(j == W)
                G(n,:) = 0;
                G(n,n) = -3;
                if(i > mid_x - (res_L/2) && i < mid_x + (res_L/2))
                    G(n,nxm) = s2;
                    G(n,nxp) = s2;
                    G(n,nym) = s2;
                    Smap(i,j) = s2;
                else
                    G(n,nxm) = s1;
                    G(n,nxp) = s1;
                    G(n,nym) = s1;
                    Smap(i,j) = s1;
                end
            else          
                G(n,:) = 0;
                G(n,n) = -4;
                %  setting my X and Y Boundaries  
                if((i > mid_x-(res_L/2) && i < mid_x+(res_L/2)) && ...
                        (j > mid_y+(res_W/2) || j < mid_y-(res_W/2)))
                
                    G(n,nxp) = s2;
                    G(n,nxm) = s2;
                    G(n,nyp) = s2;
                    G(n,nym) = s2;
                    Smap(i,j) = s2;
                else
                    G(n,nxp) = s1;
                    G(n,nxm) = s1;
                    G(n,nyp) = s1;
                    G(n,nym) = s1;
                    Smap(i,j) = s1;
                end     
            end
        end
    end
%%
    V = G\B;
    Vmap = zeros(L,W);
    for i =1:1:L
        for j =1:1:W
            n = j + (i-1)*W;
            Vmap(i,j) = V(n);
        end
    end


    E = gradient(Vmap);

    J = Smap.*E;

    area = L*W;
    J_avg = sum(sum(J))/(L*W);
    I = J_avg/area;

end