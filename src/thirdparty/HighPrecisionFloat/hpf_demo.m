%% Demo script for my High Precision Floating point decimal
% format. The user can specify any number of digits to be carried, while
% doing a variety of different numerical computations on these numbers. Not
% all MATLAB operators are defined, as my main goal here was merely to build
% a general tool for this purpose while learning to use MATLAB's OOP
% facilities. As well, it seems a useful tool to learn how one MIGHT
% accomplish the goal of a "big" float format.
%
% Considerably more information is provided in the document HPF.pdf,
% where I provide many details of the computational methods employed in
% HPF.
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com

%% Create a hpf number at the command line
% Note that you don't need to specify the number of digits in advance.
% HPF has a default of 66 total digits, 2 of which are shadowed, carried
% as guard digits, but not reported. 
DefaultNumberOfDigits 64 2
F = hpf(12)

%%
% There are 66 digits carried in this format by default, two of which are
% generally hidden from view.
F.NumberOfDigits

%%
% Note that the number 12 is an integer in MATLAB, so it is exactly
% convertable in the hpf format. Next, convert a true floating point number
% into hpf form.
F = hpf(1.2)

%%
% As we see here, the number was not truly exactly 1.2. This happens because
% I let MATLAB create the value 1.2, which is then passed into hpf and then
% converted into the decimal form for that number from the hex form carried
% internally. The digits here are the exact representation of that number as
% MATLAB has stored it in the ieee form. In fact, the last 40 digits are
% essentially random floating point trash. Also see that this
% number has a few trailing zeros.
mantissa(F)

%% 
% Had we specified F in a different way, we could have generated the exact
% decimal form. Thus
F = hpf('1.2')

%%
% This number is the exact representation of the desired decimal value.
% A nice feature of HPF is that it will parse anumber with a long string of
% decimal digits correctly. Thus we can provide a number with many digits
% and expect that hpf will store all of them exactly as we desire.
F = hpf('123.45678909876543210123456789098765432101234567890987654321012345',[100 0])

%%
% since we requested exactly 100 digits in the result, there are zeros
% appended to the end of the digit string.
mantissa(F)

%%
% Be careful though, since if you provide many digits, but then indicate
% that only a few digits be stored, then I'll do as you tell me, truncating
% the extra digits from that number.
F = hpf('123.45678909876543210123456789098765432101234567890987654321012345',[10 0])
mantissa(F)

%%
% Note that the digits of F are stored as Migits, with an initial default
% set as 4-migits
DefaultDecimalBase

%%
% We can see the migits themselves
F.Migits

%%
% Numbers are stored in a floating point form, with a separate sign and
% exponent from the Mantissa. The sign is stored as either -1, 0 or 1.
% We can extract the sign directly.
F.Sign

%%
% Extracting the exponent is also easily done.
F.Exponent

%%
% That exponent is stored as an int64 integer, so numbers as large in
% magnitude as the largest int64 number can be formed and manipulated.
intmax('int64')

%%
% HPF does include both Inf and NaN as permissable values.
hpf(inf)
hpf(nan)

%%
% In some cases a computation can have a problem. So just as 0/0 produces a NaN 
% in matlab, it does so with HPF numbers too. The latter two examples here
% create HPF NaNs.
0/0
hpf(0)/0
hpf(0)/hpf(0)

%%
% Inf is as easily generated, even by some functions. As large as the
% exponent in any floating point form can be, we can always overflow that
% value if we try hard enough.
1/0
hpf(1)/0
cot(hpf(0))
tan(hpf('pi',[64 1])/2)
exp(hpf('1e1000'))

%%
% Most numeric classes are currently supported for conversion into a hpf
% format, with the exception of complex numbers. (Sorry, perhaps in a
% future version.)

%% Special numbers in hpf
% The numbers pi and e are useful enough that I'm willing to store the
% true value for those digits, up to a realistic limit. Thus, there are
% 100000 digits available for pi, and 100000 digits stored away for e.
%
% Here I will define a 100 digit version of the number pi. Of course, I
% cannot simply type hpf(pi), as that would convert the double precision
% version of pi that MATLAB supplies into a hpf. Instead, HPF looks for
% 'pi' and 'e' as keys here.
PI = hpf('pi',100)

%% Properties of a hpf number
properties('hpf')

%% Methods currently defined for hpf numbers.
%   (This will grow with time of course.)
methods('hpf')

%% Arithmetic operations are defined on a hpf number
2*PI

%%
% So if we subtract the value of 2*pi, as defined by MATLAB, we expect a
% result that is on the order of eps.
2*PI - 2*pi

%% Many of the standard functions in mathematics have been defined for hpf numbers
% For example, I could get 100 digits in the value of e as simply
e_100 = hpf('e',100)

%%
% But I could have also found that value as
exp(hpf(1,100))

%% Changing the number of digits in an hpf number
% We can specify the number of digits to be carried in an hpf number in a
% variety of ways.
%
% First, I'll define 1/3 as a hpf number, with 100 decimal digits shown.
% in fact, there will be 104 digits stored, the last four of which are
% simply there to be conservative, for an extra bit of accuracy.
F = hpf('1',100)/3
F.Migits
F.NumberOfDigits

%%
% Reduce the number of digits in F
F = augmentdigits(F,[50,2])
F.NumberOfDigits

%%
% Of course, we can increase the number of digits stored again, but the
% information in those truncated digits is lost. Future computations will
% now be carried out with the new precision.
F = augmentdigits(F,[250,4])
F.NumberOfDigits
F.Migits

%% Trig functions
% The standard trigonometric functions are all defined (although I've
% probably missed your favorite special function as there are so many.
% I seem to be constantly adding new functions as I notice that one is
% missing.)
%
% The sine of 30 degrees = pi/6 radians is 1/2
sin(hpf('pi',200)/6)
sind(hpf('30',1000))

%%
% The cosine of 45 degrees = pi/4 radians is sqrt(2)/2
C = cosd(hpf('45',500))
% and of course, the square of that number should be 1/2
C.^2

%%
% The tangent of 45 degrees = pi/4 radians is 1
tand(hpf('pi',1000)/4)
atand(hpf('1',1000))

%% Basic arithmetic is provided, and much more.
% As long as one term in the expression is an HPF number, the computation
% is done with an HPF result. Be careful though. Here I forced 5 to be an
% hpf number, with the current DefaultNumberOfDigits, because then I took a
% square root of 5. I can add an integer to that, since HPF converts 1 on
% the fly into HPF.
DefaultNumberOfDigits 64 4
phi = (1 + sqrt(hpf(5)))/2

%%
% Beware however, that while the following result will be in HPF form, with
% the desired number of digits, that the result will be incorrect past the
% 16th digit or so, because here sqrt(5) is computed as a double precision
% number, and only then converted to an HPF form.
(hpf(1) + sqrt(5))/2

%%
% If you want to be confident that the number computed is accurate, then do
% it a second time, but with more digits of precision, using a few more
% spare digits in the number. Here, I'll subtract the more accurate result
% from phi. All the reported digits are identical to those from our
% original computation. 
phi - (1 + sqrt(hpf(5,[64 10])))/2

%% A test of some identities
% Working in 1000 decimal digits, with no shadow digits to hide any flaws...
DefaultNumberOfDigits 1000 0
x = sqrt(hpf(rand))

%%
% Should be 1
sin(x).^2 + cos(x).^2

%%
% Should be -2
(log(cosh(x)/exp(x) - 1/2) + log(hpf(2)))/x

%% Computing pi
% As I've shown, I already store the value of pi internally in HPF to a
% rather large number of digits. But, one could use HPF to compute pi
% yourself. Here I'll use Viete's formula to generate 64 digits of pi.
%   http://en.wikipedia.org/wiki/Pi
a = sqrt(hpf(2,[64 3]));
b = a/2;
tol = 10*eps(b)
niter = 0;
while 1
  niter = niter + 1;
  a = sqrt(2 + a);
  b0 = b;
  b = b*a/2;
  
  if abs(b - b0) <= tol
    break
  end
end
piest = 2/b

%%
% It took niter iterations before it had converged to my specified tolerance.
niter

%%
% As you can see, the approximation yields the same digits as those I have
% stored for pi.
hpf('pi',64)

%%
% Another choice might be Chudnovski's formula:
%   http://en.wikipedia.org/wiki/Chudnovsky_algorithm
DefaultNumberOfDigits 5000 10
DefaultDecimalBase 5
a1 = hpf(13591409); % small integers, so exact
a2 = hpf(545140134);
a3 = hpf(640320);
a4 = a3.*a3.*a3;
num = hpf(1);
den = sqrt(a4);
tol = eps(num);
piinv = hpf(0);
k = 0;
while 1
  
  piterm = num.*(a1+a2.*k)./den;
  piinv = piinv + piterm;

  if abs(piterm) < tol
    break
  end
  
  k = k + 1;
  num = num.*prod(hpf(6*k+(-5:0)));
  den = den.*a4.*(3*k-2).*(3*k-1).*(3*k).*hpf(k).^3;
  
  num = -num;

end
piest = reciprocal(12.*piinv)

%%
% The numner of iterations required for this accuracy is
k

%%
% So roughly 5000/k digits per iteration
5000/k

%%
% With an error of:
piest - hpf('pi')

%%
% Of course, we can use the rat function to give a simple rational fraction
% approximation for pi also. Perhaps you have tried rat out for doubles.
[N,D] = rat(pi)

%%
% Rat now works on HPF numbers, and as you can see, here it yields a
% rational approximation for pi accurate to 100 digits. Its not terribly
% fast, since there are a lot of divides. Sorry about that, but 100 digits
% is no problem at all.
[N,D] = rat(hpf('pi',100),hpf('1e-101'))
hpf('pi',100)
N/D

%% A timing test
% Time for a multiply in HPF should be roughly inversely quadratic
% although 35000 digits is not enough to really stress the system
DefaultNumberOfDigits 35000 0
DefaultDecimalBase 1
tic,hpf('pi')*hpf('e');toc

DefaultDecimalBase 2
tic,hpf('pi')*hpf('e');toc

DefaultDecimalBase 3
tic,hpf('pi')*hpf('e');toc

DefaultDecimalBase 4
tic,hpf('pi')*hpf('e');toc

DefaultDecimalBase 5
tic,hpf('pi')*hpf('e');toc

DefaultDecimalBase 6
tic,hpf('pi')*hpf('e');toc

%% Linear algebraic tools: chol, LU, det
% Start with a Hilbert matrix of order 5
N = 5;
k = repmat(hpf(1:N,100),N,1);
A = reciprocal(k + k' - 1);

%%
% Cholesky factor
L = chol(A);
%%
% Max error
max(max(abs(A - L'*L)))

%%
% LU factorization
[L,U,P] = lu(A);
max(max(abs(P*A - L*U)))

%%
% Determinant. For this one, lets create a little higher order matrix.
% Note that the matrix is NOT truly singular, but the deterinant is
% vanishingly tiny. Had I used too high of a order for this matrix, we
% would have seen problems even in 100 digits of precision.
N = 25;
k = repmat(hpf(1:N,100),N,1);
A = reciprocal(k + k' - 1);
det(A)

%% Use, and "Abuse" of HPF
% Of course, you can do many interesting things with the HPF toolbox. Above
% I used it to compute pi to many digits of precision. You can also use
% this tool to do computations to many digits of precision on numbers that
% you know inexactly, an abuse of both those numbers and this toolbox. As
% well, you can also use this toolbox to do computations that are better
% done using well thought out numerical methods. Good numerical analysis
% will always trump brute force computation. But, at times brute force can
% be too much of a temptation.
%
% For example, from a recent Project Euler problem (318), how many leading
% nines are there in the fractional part of (sqrt(2) + sqrt(3)).^n, where
% n is a large even number?
%
% Thus, for various values of n...
DefaultNumberOfDigits 2000 3 session
k = sqrt(hpf(2)) + sqrt(hpf(3));
for n = [2 4 8 16 32 64 128 256 512 1024]
  f = fractionalpart(k.^n)
  % How many 9's does that fractional part begin with?
  find(mantissa(f) ~= 9,1,'first') - 1
end

%%
% The problem is brute force will fail you in this quest, unless you find
% the key. That, of course, is the essence of all Project Euler problems.

