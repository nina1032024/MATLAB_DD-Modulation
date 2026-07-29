function B = get_block(Mat, blockSize, rowBlockIdx, colBlockIdx)
% =========================================================================
% Program : get_block.m
% Description :
%   Extract one square block from a block-partitioned matrix. Used for two
%   DIFFERENT partitions in Stage 2 (chap4_list.md Sim-6...Sim-9), which is
%   exactly why this takes blockSize as an explicit argument rather than
%   assuming M or N:
%
%     - G / Ȟ  are partitioned into N blocks of size M  (block index = time
%       slot n, q = m + n*M)         -> get_block(G,  M, n, n')   Sim-6, Sim-7
%     - H̃ / H  are partitioned into M blocks of size N  (block index = delay
%       m, via index n + m*N after the P permutation)
%                                    -> get_block(Htilde, N, m, m')  Sim-8, Sim-9
%
%   Block indices are 0-based (matching the book's q/m/n conventions).
%
% Input  : Mat         - square matrix, size must be an integer multiple of blockSize
%          blockSize   - block edge length
%          rowBlockIdx - 0-based row block index
%          colBlockIdx - 0-based column block index
% Output : B           - blockSize x blockSize submatrix
% =========================================================================

rows = rowBlockIdx*blockSize + (1:blockSize);
cols = colBlockIdx*blockSize + (1:blockSize);
B = Mat(rows, cols);

end
