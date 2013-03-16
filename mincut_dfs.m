function cutset = mincut_dfs(i,flow,capacity,Imx,cutset)
    cutset=[cutset;i];
    children = find(Imx(i,:));%get all children
    e=-0.1;
    for c=1:length(children)
        if sum(ismember(cutset,children(c)))>0%alread in set
            continue;
        end
        idx = Imx(i,children(c));
        if (idx<=0)
           error('not connected nodes! check it\n');
        end
        if flow(idx)-capacity(idx)<e
            cutset = mincut_dfs(children(c),flow,capacity,Imx,cutset);
        end
    end
end