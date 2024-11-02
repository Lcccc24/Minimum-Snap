function Ct = getCt(n_seg, n_order)
    %#####################################################
    % STEP 2.1: finish the expression of Ct
    %
    %
    %
    %
    %
       
    %Ct = zeros((n_order+1)*n_seg,8+4*(n_seg-1));
    %思想为列递进，找对应的行数
    for i = 1:4
        Ct(i,i)=1;
    end

    index_col = size(Ct,2);
    for i = 1:n_seg-1
        index = 8*i;
        Ct(index-4+1,index_col+i)=1;
        Ct(index+1,index_col+i)=1;
    end

    index_col = size(Ct,2);
    for i = 1:4
        index = (n_seg-1)*8+4;
        Ct(index+i,index_col+i)=1;
    end

    index_col = size(Ct,2);
    for i = 1:n_seg-1
        index = 8*i;
        index_c = index_col+(i-1)*3;
        Ct(index-4+2,index_c+1)=1;
        Ct(index+2,index_c+1)=1;
        Ct(index-4+3,index_c+2)=1;
        Ct(index+3,index_c+2)=1;
        Ct(index-4+4,index_c+3)=1;
        Ct(index+4,index_c+3)=1;
    end
end
% function Ct = getCt(n_seg, n_order)
%     % C矩阵维度为 (8*段数 )X (4+k-1+4+(k-1)*3)
%     Ct = zeros((n_order+1)*n_seg,8+(n_seg-1)*4);
%     % start,p,v,a,j，已知的起点数据
%     for i=0:3
%         Ct(i+1,i+1)=1;
%     end
%     % end,p,v,a,j，已知的终点数据
%     for i = 0:3
%         Ct(8*(n_seg-1)+4+i+1,4+(n_seg-1)+i+1)=1;
%     end
%     if n_seg~=1
%         for i=0:n_seg-2
%             % p，已知的位置
%             Ct(8*i+4+1,4+i+1)=1;
%             Ct(8*(i+1)+1,4+i+1)=1;
%             % v，未知的速度
%             Ct(8*i+4+2,8+(n_seg-1)+i*3+1)=1;
%             Ct(8*(i+1)+2,8+(n_seg-1)+i*3+1)=1;
%             % a，未知的加速度
%             Ct(8*i+4+3,8+(n_seg-1)+i*3+2)=1;
%             Ct(8*(i+1)+3,8+(n_seg-1)+i*3+2)=1;
%             % j，未知的jerk
%             Ct(8*i+4+4,8+(n_seg-1)+i*3+3)=1;
%             Ct(8*(i+1)+4,8+(n_seg-1)+i*3+3)=1;
%         end
%     end
% end
