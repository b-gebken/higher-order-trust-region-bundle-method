% Adds elements to the memory-struct. If the maximum size of the memory is
% reached, then the oldest elements are removed.  

function memory = add_to_memory(new_sample_pts,new_oracle_vals,memory)

    q = size(new_oracle_vals,1)-1;

    memory.sample_pts = [memory.sample_pts,new_sample_pts];
    memory.oracle_vals = [memory.oracle_vals,new_oracle_vals];

    if(size(memory.sample_pts,2) > memory.max_size)
        memory.sample_pts = memory.sample_pts(:,end-memory.max_size+1:end);
        memory.oracle_vals = memory.oracle_vals(:,end-memory.max_size+1:end);
    end
end

