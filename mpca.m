function [D, V] = mpca(X) 
  m = mean(X);
  Xc = bsxfun(@minus, X ,m);
  E = cov(Xc);
  [V,D] = eig(E);
  [D, i] = sort(diag(D), 'descend');
  V = V(:,i);
end