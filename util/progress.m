function progress(cur_iter, num_iters)
persistent old_tmp;

if cur_iter == 1
    msg = sprintf('%.3f%% done.',cur_iter/num_iters*100);
    old_tmp = msg;
elseif cur_iter == num_iters
    tmp = sprintf('%.3f%% done.',cur_iter/num_iters*100);
    msg =  [char(8)*ones(1,length(old_tmp)),tmp];
    old_tmp = tmp;
    fprintf(['%' num2str(length(msg)) 's\n'],msg);
    return;
else
    tmp = sprintf('%.3f%% done.',cur_iter/num_iters*100);
    msg =  [char(8)*ones(1,length(old_tmp)),tmp];
    old_tmp = tmp;
end
fprintf(['%' num2str(length(msg)) 's'],msg);

end