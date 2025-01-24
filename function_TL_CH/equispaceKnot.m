function [xsi_vec_new] = equispaceKnot(xsi_vec,num_points)

% Ordina e rimuove i duplicati
unique_xsi = unique(xsi_vec);

% Controlla se il vettore iniziale contiene solo 0 e 1
    if isequal(unique_xsi, [0 1])
        % Crea punti equidistanti tra 0 e 1 usando linspace
        xsi_vec_new = linspace(1/(num_points+1), 1-1/(num_points+1), num_points);
    else
        % Inizializza i nuovi punti
        new_points = [];
    for p = 1:num_points
            % Calcola i segmenti tra i punti
            segments = diff(unique_xsi);
    
            % Trova il segmento più lungo
            [~, idx_max] = max(segments);
    
            % Aggiungi un punto al centro del segmento più lungo
            new_point = unique_xsi(idx_max) + segments(idx_max) / 2;
            
            % Inserisci il nuovo punto nella lista
            unique_xsi = [unique_xsi, new_point];
            
            % Rimuovi duplicati
            unique_xsi = unique(unique_xsi);
            
            % Salva il nuovo punto
            new_points = [new_points, new_point];
    end
        
        % Ordina tutti i punti
        xsi_vec_new = sort(new_points);
    end
end