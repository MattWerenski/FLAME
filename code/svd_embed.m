function x = svd_embed(walks, ndim)
  [nnetworks, ngene, ~]  = size(walks);
  RR_sum = zeros(ngene);
  for i = 1:nnetworks
    W = walks(i);
    R = log(W + 1/ngene); % smoothing
    RR_sum = RR_sum + R * R';
  end
  clear R Q A
  [V, d] = eigs(RR_sum, ndim);
  x = diag(sqrt(sqrt(diag(d)))) * V';
end