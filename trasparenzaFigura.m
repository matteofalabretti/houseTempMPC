function trasparenzaFigura(trasparenza)
    % Supponiamo che la figura esistente sia la figura corrente
    fig = gcf;
    
    % Trova tutti gli oggetti Patch nella figura
    patches = findall(fig, 'Type', 'Patch');
    
    % Modifica la trasparenza di tutti i Patch trovati
    for k = 1:length(patches)
        set(patches(k), 'FaceAlpha', trasparenza); % Imposta la trasparenza a 0.1
    end
    
end

