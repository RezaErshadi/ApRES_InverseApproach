function misfit = FUNC_GetTheMisfit(ObsC,EstParC,CostFunc,nrm)
    a = ObsC{CostFunc}(:);
    b = EstParC{CostFunc}(:);
    a(isinf(a)) = nan;
    b(isinf(b)) = nan;
%     if CF == 16
%         a = unwrap(a);
%         b = unwrap(b);
%     end
    switch nrm
        case 0
            misfit = a - b;
        case 1
            nrm_a = (a-nanmean(a(:))) ./ (nanstd(a(:)));
            nrm_b = (b-nanmean(b(:))) ./ (nanstd(b(:)));
            misfit = nrm_a - nrm_b;
    end 
    misfit(isnan(misfit)) = [];
end