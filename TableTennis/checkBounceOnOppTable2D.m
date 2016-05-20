function bounce = checkBounceOnOppTable2D(pos,table)

dist_to_table = table.DIST;
table_length = table.LENGTH;
table_z = table.Z;

bounce = false;
oppTableCenterY = dist_to_table - 3*table_length/4;

% if filter state goes below threshold bounce is highly likely
tol = 2e-2;
if pos(2) < table_z + tol && abs(pos(1) - oppTableCenterY) < table_length/4
    % if it bounces do not hit
    disp('Ball bounces on opp table! Not hitting!');
    bounce = true;
end