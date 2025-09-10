% Author: Ettore Napoli, University of Salerno
% Date:   2025-08
%
% Test alphaevolve algorithm on random matrices
% tol     : tolerance onthe difference of each computed element.
% 		    Since the algorithm in this phase uses real numbers there are small
% 		    differences in the results that are not errors but due to numerical precision.
% n_tests : The number of random matrices multiplied. 
% 			Change the value to more thotoughly test the algorithm
%
% -------------------------------------------------------------------------
% License:
% This code is released under the MIT License.
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% -------------------------------------------------------------------------

tol = 1e-3;
n_tests = 1000;
error_count = 0;

for k = 1:n_tests
    A = randn(4) + 1i * randn(4);
    B = randn(4) + 1i * randn(4);
    
    C1 = A * B;
    C2 = alphaevolve_4x4_complex(A, B);
    
    diff_C1_C2 = abs(C1 - C2);
    
    if any(diff_C1_C2(:) > tol)
        fprintf('Test %d: Error, C1 and C2 differ (max diff = %.3e)\n', k, max(diff_C1_C2(:)));
        error_count = error_count + 1;
    end
end

if error_count == 0
    disp('All tests passed within tolerance.');
else
    fprintf('%d tests failed with errors > %.1e.\n', error_count, tol);
end
