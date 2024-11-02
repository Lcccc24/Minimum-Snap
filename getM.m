
function M = getM(n_seg, n_order, ts)
    M = [];
    for k = 1:n_seg
        M_k = zeros(8, n_order+1);
        %#####################################################
        % STEP 1.1: calculate M_k of the k-th segment 
        for r=0:3
            for c=r:n_order
                M_k(r+1,c+1)=(factorial(c)/factorial(c-r))*0^(c-r);
                M_k(r+1+4,c+1)=(factorial(c)/factorial(c-r))*ts(k)^(c-r);
            end
        end
        M = blkdiag(M, M_k);
    end
end
