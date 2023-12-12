function di = DistanceROI(Node, ROI) % Определение массива расстояний между точкой и массивом ROI
    di = sqrt((Node(1)-ROI(:,1)).^2+(Node(2)-ROI(:,2)).^2+(Node(3)-ROI(:,3)).^2);
end