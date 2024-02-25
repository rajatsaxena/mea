function channelmap = AP_loadChannelNum(probemap, channel_positions)
    channel_map = [probemap.xcoords probemap.ycoords];
    % Initialize array to store indices
    matching_indices = [];

    % Iterate through each row of array1
    for i = 1:size(channel_positions, 1)
        % Compare each row of array1 with all rows of array2
        % Check if the rows are identical using the ismember function
        [~, idx] = ismember(channel_positions(i, :), channel_map, 'rows');

        % Find the indices where the row of array1 matches array2
        matching_indices = [matching_indices; idx];
    end
    
    % return matching channel numbers
    chnum = probemap.chanMap(matching_indices);
    shnum = zeros([length(chnum) 1]);
    shnum(channel_positions(:,1)<=100,1) = 1;
    shnum(channel_positions(:,1)>100,1) = 2;
    channelmap=zeros([length(chnum) 4]);
    channelmap(:,1)=chnum; % channel number
    channelmap(:,2)=shnum;
    channelmap(:,3)=channel_positions(:,1); % x coords
    channelmap(:,4)=channel_positions(:,2); % y coords

end