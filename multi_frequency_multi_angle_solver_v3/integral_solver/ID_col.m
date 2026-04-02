function [Js,Z,k] = ID_col(A,acc)
% This program computes row ID with a given accuracy acc
    normA = norm(A,'f');
    [m,n] = size(A);
    dim = min(m,n);
    
    if dim <= 2000
        ss = svd(A);
        s = 1;
        k  = max(sum(ss/ss(1) > .1*acc*normA),s);
        [~, R, ind] = qr(A,0);
        T = R(1:k,1:k) \ R(1:k, (k+1):size(R,2));

        Z = zeros(k,n); Z(:,ind) = [eye(k), T];
        Js = ind(1:k);
    end

    if dim > 2000 % use randomized method
        disp('randomized method is used.')
        [~,B,ss] = random_qb(A,acc);
        s = 1;
        k  = max(sum(ss/ss(1) > .1*acc*norm(A)),s);
        [~, R, ind] = qr(B,0);
        T = R(1:k,1:k) \ R(1:k, (k+1):size(R,2));

        Z = zeros(k,n); Z(:,ind) = [eye(k), T];
        Js = ind(1:k);
    end

end

function [Q,B,ss] = random_qb(A,acc)
    normA = norm(A,'f');
    [m,n] = size(A);
    s = min(m,n);
    b = min(450,s);
    tol = 0.1*acc*normA*sqrt(s); % a generous indication to see if exit

    Q = []; B = [];
    for j = 1:s
        G = randn(n,b);
        [U,~] = qr(A*G,0);
       if j >1
         [U,~] = qr(U - Q*(Q'*U),0); % orthogonalize new U to existing Q
       end
       B = [B; U'*A];
       Q = [Q, U];
       A = A - U*B(end-b+1:end,:); % substract the newest U*B
       %norm(A,'fro')

        if norm(A,'fro') < tol
            ss = svd(B);
            k = round(.97*length(ss)); 
            if ss(k) > 0.1*normA*acc  % make sure the singular value satisfies the acc
              tol = tol/2;
            else
              break
            end
        end
    end

end
