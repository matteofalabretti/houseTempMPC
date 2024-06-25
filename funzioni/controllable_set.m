function [H_nsteps, h_nsteps , Np] = controllable_set(Hx, hx, Hu, hu, G, g, A, B, point)

    % Calcolo del n-step controllable_set di un sistema lineare
%   Input:
%       - Hx*x <= hx: vincoli sullo stato
%       - Hu*u <= hu: vincoli sull'ingresso
%       - G <= g: la regione obbiettivo
%       - A, B: matrici del sistema
%       - point: il punto du partenza del nostro sistema

    n = size(A,2);
    
    % candidato iniziale
    H_ii_steps = G;
    h_ii_steps = g;
    tic;
    for i = 1:1000

        % Calcoliamo il set ad un passo rispetto a quello precedente
        % x_next = Ax + Bu <= [hx;hu]
        temp = Polyhedron('A', [H_ii_steps*A, H_ii_steps*B; zeros(size(Hu, 1), n), Hu], 'b', [h_ii_steps;hu]);
    
        % Proiezioni in R^n    
        % temp = temp.computeVRep();
        temp = projection(temp, 1:n);
        % temp = temp.minHRep();
        % temp = temp.computeHRep();

        % Intersezione con X := {x | Hx*x <= hx}
        H_ii_steps = [temp.A;Hx];
        h_ii_steps = [temp.b;hx];

        
        disp("Fine Calcolo Passo " + i )
        
        % Verifichiamo se il Poliedro contiene il punto di partenza
        temp = Polyhedron( 'A' , H_ii_steps , 'b' ,h_ii_steps);
    
        if temp.contains(point)
            % se contiene il punto usciamo dal ciclo
            break
        end
    
    end
    
    disp("Tempo impiegato: " + toc);

    temp = temp.minHRep();
    Np = i;
    H_nsteps = temp.A;
    h_nsteps = temp.b;
    
    
end