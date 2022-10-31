function lists = data_association_listoverhaul(lists,z_ID,next_idx,Rnext_idx,params)

% z_ID is a column of the detection matrix

% The feature ID should correspond to the feature's index in the cell array
% features{ID}.status should be 0 for uninitialized features, 1 for
%   previously-initialized features (potential for more later)

% In the current detection column, which features have a value of 1 instead
% of zero? These are the features that are currently visible.
detections = find(z_ID > 0);
n_detect = length(detections);
deactivate_list = [];
active_list = lists.active;
n_parent = lists.n_parent;
lists.new = [];

if ~isempty(active_list)

    % Figure out which features are already active and were re-detected
    [~,~,ib] = intersect(detections,active_list(:,1),'stable');
    persist_list = active_list(ib,:);

    if size(persist_list,1) < params.n_max

        new_candidates = setdiff(detections,persist_list(:,1));
        [~,id] = setdiff(active_list(:,1),detections);
        deactivate_list = active_list(id,:);
        update_list = persist_list;

        if ~isempty(new_candidates)

            j = 1;
            psize = size(persist_list,1);
            while psize < params.n_max

                next_feature.ID = new_candidates(j,1);
                next_feature.idx = next_idx;
                next_feature.Ridx = Rnext_idx;
                next_feature.parent = n_parent + 1;
                next_feature.status = 0;
                
                lists.features{next_feature.ID} = next_feature;
                update_list = [update_list;next_feature.ID];
                lists.new = [lists.new;next_feature.ID];
                
                next_idx = next_idx + 1;
                Rnext_idx = Rnext_idx + 1;
                psize = psize + 1;
                j = j + 1;

                if j > size(new_candidates,1)
                    break
                end

            end

        end

    else

        update_list = persist_list;

    end

else

    % Activate some features
    if n_detect < params.n_max

        update_list = detections;

    elseif n_detect >= params.n_max

        update_list = detections(1:params.n_max);

    end

    for j = 1:size(update_list,1)

        next_feature.ID = update_list(j,1);
        next_feature.idx = next_idx;
        next_feature.Ridx = Rnext_idx;
        next_feature.parent = n_parent + 1;
        next_feature.status = 0;

        lists.features{next_feature.ID} = next_feature;
        lists.new = [lists.new;next_feature.ID];

        next_idx = next_idx + 1;
        Rnext_idx = Rnext_idx + 1;

    end

end

lists.active = update_list;
lists.deactivate = deactivate_list;

end