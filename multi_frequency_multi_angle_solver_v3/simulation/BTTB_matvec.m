function y = BTTB_matvec(Cpp, x) 
% Fast calculation of convolution by FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% Cpp: Green's function convolution operator without k^2 prefactor (k^2 is
% contained in P), size (2*ny-1, 2*nx-1).
% x: polarization currents, with k^2 prefactor, size (ny*nx, 1)

% Outputs:
% y: G*x, convolution results, size (ny*nx,1)

% Note:
% Escat = G*P, G is Green's function convolution matrix, P is polarization
% currents, size (ny*nx,1), G is a Toeplitz matrix with Toeplitz Block (BTTB), 
% size (ny*nx, ny*nx). This direct convolution is slow. There is a fast algorithm to calculate
% the matrix vector multiplication for a BTTB matrix
% We first construct a Block circulant matrix with circulant blocks (BCCB) C,
% (2*nx-1,2*nx-1) blocks, each block is a (2*ny-1,2*ny-1) matrix.
% G*P = C*P_padded, P_padded add zero paddings to P
% BCCB matrix vector multiplication can be computed by FFT
% A BCCB matrix is completely determined by its first column, which has
% (2*ny-1)*(2*nx-1) elements, we can reshape them into a (2*ny-1, 2*nx-1)
% matrix. This is the matrix Cpp that we get from build_greens_function
% For more details, you can ask Zhaowei or Zeyu for notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Escat = G*P, G is a BTTB, Block Toeplitz matrix with Toeplitz Block
[N, M] = size(Cpp); % dimensions
n = (N+1)/2;
m = (M+1)/2;

% convert the nm by 1 input vector into a n by m matrix and add zero
% paddings to stuff it to 2n-1 by 2m-1
X = reshape(x, n, m);
Xpp = zeros(N, M); % X with zero paddings
Xpp(1:n, 1:m) = X;

% get output and strip down the paddings, reshape the resulting matrix 
% to vector form 
Ypp = ifft2(Cpp .* fft2(Xpp));
Y = Ypp(1:n, 1:m);
y = reshape(Y, n*m, 1);
end