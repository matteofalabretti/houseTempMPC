function [A_cal , A_cal_n ,B_cal , B_cal_n , Q_cal , R_cal] = Calligrafica(A , B , Q , R , S , N)
    % Funzione che calcola le matrici nella forma raccolta
    
    % Inizzializzazioni matrici
    A_cal = [];
    B_cal = [];
    
    for i = 1:N
        mat = A^i;
        % Calcolo di A_cal
        A_cal = [A_cal;mat];
        if i == N
            A_cal_n = mat;
        end
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
    
        if i == N
            B_cal_n = riga_B_cal;
        end

    end
    
    %Calcolo di Q_cal
    Q_cal = diagonaleEccettoUltima(Q , S , N);
    
    %Calcolo di R cal
    R_cal = kron(eye(N) , R );
    
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