function subTotal = computeInputEffect(A_d, B_d, mu, p, time)
   n = length(A_d);
   subTotal = zeros(n, p);
   for s=1:time
    subTotal = subTotal + A_d^(s-1)*B_d*mu(:, (time-s)*p+1:(time-s+1)*p);
   end
end