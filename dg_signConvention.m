function [U,V] = dg_signConvention(Lz,Rz,U,V)

for mIdx = 1:size(V,2)
    tL = Lz * U(:,mIdx);
    [~, domIdx] = max(abs(U(:,mIdx)));
    L_dom = Lz(:, domIdx);
    if corr(tL, L_dom) < 0
        U(:,mIdx) = -U(:,mIdx);
        V(:,mIdx) = -V(:,mIdx);
    end
end
