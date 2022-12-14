function outsig = overlap_add(frames,idx_mat)
    [num_frames] = size(idx_mat);
    overadd_count = zeros(1, max(idx_mat(:)));
    outsig = zeros(1, max(idx_mat(:)));

    for n = 1:num_frames
        outsig(idx_mat(n,:)) = outsig(idx_mat(n,:)) + frames(n,:);
        overadd_count(idx_mat(n,:)) = overadd_count(idx_mat(n,:)) + 1;
    end

    % averaging the overlapped estimates for each sample
    outsig = outsig./overadd_count;
end