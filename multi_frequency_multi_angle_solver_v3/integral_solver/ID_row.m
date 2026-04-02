function [Js,X,out] = ID_row(A,acc)
% This program computes row ID with a given accuracy acc
% A = X*R, A is m by n, X is m by k, R is k by n
% X(Js,:) = eye(k), so Js tells identical rows of A in X
% Because of the above reason, R needs not to be computed
%
                % A \approx Z*A(Js,:)
                % Z(Js,:) is k by k idenitty matrix
% out is an optional output to indicate the compressed rank k
    
    if nargout == 3
        [Js,X,k] = ID_col(A.',acc);
        X = X.';
        out = k;
    elseif nargout == 2
        [Js,X] = ID_col(A.',acc);
        X = X.';
    end
end

