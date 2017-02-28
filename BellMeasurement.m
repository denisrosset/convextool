function BM = BellMeasurement
% BellMeasurement Returns a 1x4 cell array containing 4x4 Herm. matrices
    BM = cell(1, 4);
    BM{1} = [1 0 0 1 % phi plus
             0 0 0 0
             0 0 0 0
             1 0 0 1]/2;
    BM{2} = [ 1 0 0 -1 % phi minus
              0 0 0  0
              0 0 0  0
             -1 0 0  1]/2;
    BM{3} = [0 0 0 0
             0 1 1 0
             0 1 1 0
             0 0 0 0]/2;
    BM{4} = [0  0  0 0
             0  1 -1 0
             0 -1  1 0
             0  0  0 0]/2;
end
