function X=bandf(X, Fc1, Fc2)
    X = X - mean(X');
    [B] = fir1(32, Fc1, 'low');
    XX=filtfilt(B, 1, X);
    X = resample(XX, 1, 2);
    X = resample(X, 1, 2);
    [B1] = fir1(8, Fc2*4, 'low');
    XXX = filtfilt(B1, 1, X);
    X = resample(XXX, 2, 1);
    X = resample(X, 2, 1);
    [n,N] = size(XX);
    X=XX-X(:,1:N);
end
