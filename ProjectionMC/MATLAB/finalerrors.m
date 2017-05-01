function [E0_final, E0error_final, x0_final, x0error_final] = finalerrors(E0, E0error, x, bx, pinterval)



%% Calcs
%Determine prange
prange = [find(pinterval>=110, 1, 'last'):find(pinterval<=90, 1)];

%Determine the final energy and its error
E0_final = mean(E0(prange));
E0error_final = sqrt(dot(E0error(prange),E0error(prange))/length(prange));

%Determine the final position and its error
x0_vec = zeros(1,21);
x0error_vec = zeros(1,21);
index = 1;
for p = 90:110
    [x0_vec(index), x0error_vec(index)] = meanposition( x, bx, p+1 );
    index = index+1;
end
x0_final = mean(x0_vec);
x0error_final = sqrt(dot(x0error_vec,x0error_vec)/length(prange));



end