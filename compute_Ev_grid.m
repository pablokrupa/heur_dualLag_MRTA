%% This function computes the adjacency matrix Ev of an nxm grid

function Ev = compute_Ev_grid(n, m)
    dim = n*m;
    Ev = zeros(dim, dim);
    for r = 0:(n-1)
        for c = 0:(m-1)
            i = r*m + c;
            if c > 0
                Ev(i,i+1) = 1;
                Ev(i+1,i) = 1;
            end
            if r > 0
                Ev(i-m+1,i+1) = 1;
                Ev(i+1,i-m+1) = 1;
            end
        end
    end
end