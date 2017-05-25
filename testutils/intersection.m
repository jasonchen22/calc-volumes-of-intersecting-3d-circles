function p = intersection(cc1, nn1, cc2, nn2, cc3, nn3)

D_(1,:) = nn1; D_(2,:) = nn2; D_(3,:) = nn3;
dd = det(D_);
if dd == 0, p = []; return; end

b1 = dot(nn1, cc1); b2 = dot(nn2, cc2); b3 = dot(nn3, cc3);
p = D_\[b1 b2 b3]';
