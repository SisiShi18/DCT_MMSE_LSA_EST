function delta = time_diff(coef, theta)
    % c: coefficients(ie.,MFCC)
    % theta: delta window
    [num_fra, num_coef] = size(coef); % no. of frames, and dimension of the coefficient.
    vf = (theta:-1:-theta);
    vf = vf/sum(vf.^2);
    ww = ones(theta, 1);
    cx = [coef(ww,:); coef; coef(num_fra*ww,:)];
    vx = reshape(filter(vf, 1, cx(:)), num_fra+2*theta, num_coef);
    vx(1:2*theta,:) = [];
    delta = vx;
end