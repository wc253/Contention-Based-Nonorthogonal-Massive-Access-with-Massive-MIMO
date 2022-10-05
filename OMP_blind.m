function [B, H] = OMP_blind(Y, P,SNR)
a=1.2;%阈值调整参数
max_itr_OMP = max(size(P));
B = [];
Y_r = Y;
t1 = 1;
[~,n]=size(P);
H = zeros(n,size(Y,2));
norm_save(t1) = norm(Y_r,'fro');
while 1
    t1 = t1 + 1;
    B_last = B;    
    [~,k] = max(sum(abs(P'*Y_r).^2,2));
    B = union(k,B);
    H(B,:) = P(:,B) \ Y;
    Y_r = Y - P(:,B)*H(B,:);
    norm_save(t1) = norm(Y_r,'fro');
    
    if 20*log10( norm(P(:,B)*H(B,:),'fro')/norm_save(t1)) > a*SNR
        break;
    end    
       
    if norm_save(t1)/norm_save(t1-1) >= 1 
        B = B_last';
        H(B,:) = pinv(P(:,B))*Y;
        break;
    end
    
    if t1 >max_itr_OMP
        break;
    end
    
end
