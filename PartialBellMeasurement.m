function BM = PartialBellMeasurement
% BellMeasurement Returns a 1x4 cell array containing 4x4 Herm. matrices
    BM = cell(1, 2);
    BM{1} = [1 0 0 1 % phi plus
             0 0 0 0
             0 0 0 0
             1 0 0 1]/2;
    BM{2} = eye(4) - BM{1};
end
