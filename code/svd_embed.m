function x = svd_embed(walks, ndim)
  ngene = shape(walk, 2);
  RR_sum = zeros(ngene);
  for i = 1:length(network_files)
    W = walks(i);
    R = log(W + 1/ngene); % smoothing
    RR_sum = RR_sum + R * R';
  end
  clear R Q A
  [V, d] = eigs(RR_sum, ndim);
  x = diag(sqrt(sqrt(diag(d)))) * V';
end