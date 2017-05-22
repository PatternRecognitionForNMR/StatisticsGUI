function [D, W_lda] = lda(X,y)
classes = unique(y);
% prealloc for sum
Sw = zeros(size(X,2)); 
Sb = zeros(size(X,2));
mu = mean(X);

for i = 1:length(classes)
    Xi = X(find(y == classes(i)),:);
    n = size(Xi,1);
    mu_i = mean(Xi);
    XMi = bsxfun(@minus, Xi, mu_i);
    Sw = Sw + (XMi' * XMi );
    MiM =  mu_i - mu;
    Sb = Sb + n * MiM' * MiM;
end

[W_lda, D] = eig(Sw\Sb); % looks wrong but is correct, try inv(Sw)*Sb if in doubt
[D, i] = sort(diag(D), 'descend');
W_lda = real(W_lda(:,i));
end

% function [D, W_lda] = lda(X,y)
% dimension = size(X,2);
% labels = unique(y);
% C = length(labels);
% Sw = zeros(dimension,dimension);
% Sb = zeros(dimension,dimension);
% mu = mean(X);
% 
% for i = 1:C
%     Xi = X(find(y == labels(i)),:);
%     n = size(Xi,1);
%     mu_i = mean(Xi);
%     XMi = bsxfun(@minus, Xi, mu_i);
%     Sw = Sw + (XMi' * XMi );
%     MiM =  mu_i - mu;
%     Sb = Sb + n * MiM' * MiM;
% end
% 
% [W_lda, D] = eig(Sw\Sb);
% [D, i] = sort(diag(D), 'descend');
% W_lda = real(W_lda(:,i));
% end
