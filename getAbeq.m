function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % p,v,a,j constraint in start, 
    Aeq_start = zeros(4, n_all_poly);
    beq_start = zeros(4, 1);
    % STEP 2.1: write expression of Aeq_start and beq_start
    %
    for k = 0:3
        Aeq_start(k+1,k+1) = factorial(k);
    end
    beq_start = start_cond';
    %
    %
    
    %#####################################################
    % p,v,a constraint in end
    Aeq_end = zeros(4, n_all_poly);
    beq_end = zeros(4, 1);
    % STEP 2.2: write expression of Aeq_end and beq_end
    %
    for k = 0:3
        for j = k:n_order
            Aeq_end(k+1,(n_order+1)*(n_seg-1)+j+1) = factorial(j)/factorial(j-k)*(ts(n_seg)^(j-k));
        end
    end
    beq_end = end_cond';
    %
    %
    
    %#####################################################
    % position constrain in all middle waypoints
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    % STEP 2.3: write expression of Aeq_wp and beq_wp
    %

    for k = 0:n_seg-2
        Aeq_wp(k+1,(n_order+1)*(k+1)+1) = 1;
        beq_wp(k+1,1) = waypoints(k+2);
    end
    %
    %
    
    %#####################################################
    % position continuity constrain between each 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);
    % STEP 2.4: write expression of Aeq_con_p and beq_con_p
    %
    index= 0; %表示位置
    for k = 0:n_seg-2
        for j = index:n_order
            Aeq_con_p(k+1,k*(n_order+1)+j+1)=factorial(j)/factorial(j-index)*(ts(k+1)^(j-index));
        end

        for j = index:n_order
            Aeq_con_p(k+1,(k+1)*(n_order+1)+j+1)=-factorial(j)/factorial(j-index)*(0^(j-index));
        end
    end

            

    
    %#####################################################
    % velocity continuity constrain between each 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);
    % STEP 2.5: write expression of Aeq_con_v and beq_con_v
    %
    %
    %
    %
    index= 1; %表示速度
    for k = 0:n_seg-2
        for j = index:n_order
            Aeq_con_v(k+1,k*(n_order+1)+j+1)=factorial(j)/factorial(j-index)*(ts(k+1)^(j-index));
        end

        for j = index:n_order
            Aeq_con_v(k+1,(k+1)*(n_order+1)+j+1)=-factorial(j)/factorial(j-index)*(0^(j-index));
        end
    end
    %#####################################################
    % acceleration continuity constrain between each 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);
    % STEP 2.6: write expression of Aeq_con_a and beq_con_a
    %
    %
    %
    %

    index= 2; %表示加速度
    for k = 0:n_seg-2
        for j = index:n_order
            Aeq_con_a(k+1,k*(n_order+1)+j+1)=factorial(j)/factorial(j-index)*(ts(k+1)^(j-index));
        end

        for j = index:n_order
            Aeq_con_a(k+1,(k+1)*(n_order+1)+j+1)=-factorial(j)/factorial(j-index)*(0^(j-index));
        end
    end
    
    %#####################################################
    % jerk continuity constrain between each 2 segments
    Aeq_con_j = zeros(n_seg-1, n_all_poly);
    beq_con_j = zeros(n_seg-1, 1);
    % STEP 2.7: write expression of Aeq_con_j and beq_con_j
    %
    %
    %
    %
    index= 3; %表示jerk
    for k = 0:n_seg-2
        for j = index:n_order
            Aeq_con_j(k+1,k*(n_order+1)+j+1)=factorial(j)/factorial(j-index)*(ts(k+1)^(j-index));
        end

        for j = index:n_order
            Aeq_con_j(k+1,(k+1)*(n_order+1)+j+1)=-factorial(j)/factorial(j-index)*(0^(j-index));
        end
    end
    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
    beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
end

% function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
%     n_all_poly = n_seg*(n_order+1);
%     %#####################################################
%     % p,v,a,j constraint in start, 
%     Aeq_start = zeros(4, n_all_poly);
%     beq_start = zeros(4, 1);
%     % STEP 2.1: write expression of Aeq_start and beq_start
%     for k= 0:3
%         Aeq_start(k+1,k+1) = factorial(k);
%     end
%     beq_start = start_cond'; 
%     %#####################################################
%     % p,v,a constraint in end
%     Aeq_end = zeros(4, n_all_poly);
%     beq_end = zeros(4, 1);
%     % STEP 2.2: write expression of Aeq_end and beq_end
%     start_idx_2 = (n_order + 1)*(n_seg - 1);
%     for k=0 : 3
%         for i=k : 7
%             Aeq_end(k+1,start_idx_2 + 1 + i ) = factorial(i)/factorial(i-k)*ts(n_seg)^(i-k);
%         end
%     end
%     beq_end = end_cond';
%     
%     %#####################################################
%     % position constrain in all middle waypoints
%     Aeq_wp = zeros(n_seg-1,n_all_poly);
%     beq_wp = zeros(n_seg-1,1);
%     for i=0:n_seg-2
%         start_idx_2 = (n_order + 1)*(i+1);
%         Aeq_wp(i+1,start_idx_2+1) = 1;
%         beq_wp(i+1,1) = waypoints(i+2);
%     end
%     
%     Aeq_con = zeros((n_seg-1)*4, n_all_poly);
%     beq_con = zeros((n_seg-1)*4, 1);
%     for k=0:3
%         for j=0:n_seg-2
%             for i = k:7
%                 start_idx_1 = (n_seg-1)*k;
%                 start_idx_2 = (n_order+1)*j;
%                 Aeq_con(start_idx_1 + j + 1,start_idx_2 + i+1)=...
%                     factorial(i)/factorial(i-k)*ts(j+1)^(i-k);
%                 if(i == k)
%                     Aeq_con(start_idx_1+j+1,start_idx_2+(n_order+1)+i+1) = ...
%                                                             -factorial(i);
%                 end
%             end
%         end
%     end
%     
%     
%     %#####################################################
%     % combine all components to form Aeq and beq   
%     Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
%     beq = [beq_start; beq_end; beq_wp; beq_con];
% end