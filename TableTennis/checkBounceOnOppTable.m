function bounce = checkBounceOnOppTable(filter,table)

dist_to_table = table.DIST;
table_length = table.LENGTH;
table_z = table.Z;
table_width = table.WIDTH;

bounce = false;
oppTableCenterY = dist_to_table - 3*table_length/4;

% if filter state goes below threshold bounce is highly likely
tol = 2e-2;
if filter.x(3) < table_z + tol && abs(filter.x(1)) < table_width/2 && ...
    abs(filter.x(2) - oppTableCenterY) < table_length/4
    % if it bounces do not hit
    disp('Ball bounces on opp table! Not hitting!');
    bounce = true;
end