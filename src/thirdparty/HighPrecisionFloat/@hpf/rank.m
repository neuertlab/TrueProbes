function r = rank(A,tol)
%RANK   Matrix rank, for HPF arrays
%   RANK(A) provides an estimate of the number of linearly
%   independent rows or columns of a matrix A.
%
%   RANK(A,TOL) is obtained using rref, as modified for HPF.
%
%   Class support for input A:
%      float: double, single

%   Copyright 1984-2015 The MathWorks, Inc.

% modfiied for HPF to use rref
if nargin == 1
  [~,jb] = rref(A);
else
  [~,jb] = rref(A,tol);
end

% the rank is just the number of pivot columns found by rref
% using the default tolerance, or using tol as supplied by
% the user here.
r = numel(jb);

