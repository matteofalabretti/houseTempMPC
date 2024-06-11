function [A_cal , B_cal , Q_cal , R_cal , G_cal , g_cal,  Hx_cal , hx_cal , Hu_cal , hu_cal] = Calligrafica(A , B , Q , R , S , N , Hx , hx , Hu , hu , G , g)
    % Funzione che calcola le matrici nella forma raccolta
    
    % Inizzializzazioni matrici
    A_cal = [];
    B_cal = [];
    
    for i = 1:N
        % Calcolo di A_cal
        A_cal = [A_cal;A^i];
    end
    
    % Calcolo di B cal
    for i = 1:N
        riga_B_cal = [];
        for j = 1:N
            if (i-j) < 0
                riga_B_cal = [riga_B_cal , zeros(height(B) , width(B))];
            else
                riga_B_cal= [riga_B_cal , A^(i-j) * B];
            end
        end
    
        B_cal = [B_cal ; riga_B_cal];
    end
    
    %Calcolo di Q_cal
    Q_cal = diagonaleEccettoUltima(Q , S , N);
    
    %Calcolo di R cal
    R_cal = kron(eye(N) , R );
    
    %calcolo di G cal
    matrici = {};
    for i = i:N
        if(i = N-1)
            matrici{i} = G;
        else
            matrici{i} = eye();
        end
    end
    G_cal = blkdiag(matrici(:));

    %% Non serve
    % %Cacloclo di Hx_cal
    % Hx_cal = diagonaleEccettoUltima(Hx , G , N);
    % 
    % %Cacloclo di hx_cal
    % hx_cal = [];
    % for i = 1:N-1
    %    hx_cal = [hx_cal;hx];
    % end
    % hx_cal = [hx_cal;g];
    % 
    % % Calcolo di Hu_cal
    % Hu_cal = kron(eye(N-1) , Hu);
    % 
    % %calcolo di hu_cal
    % hu_cal = [];
    % for i = 1:N-1
    %    hu_cal = [hu_cal;hu];
    % end
end

%Funzione che crea una matrice diagonale a blocchi con tutti i blocchi
%uguali a p eccetto per l'ultima che è s, con N il numero di blocchi
function mat = diagonaleEccettoUltima(p , s , N)
    
    matrici = cell(1 , N);
    
    for i = 1:(N-1)
        matrici{i} = p;
    end
    matrici{N} = s;
    mat = blkdiag(matrici{:});
end