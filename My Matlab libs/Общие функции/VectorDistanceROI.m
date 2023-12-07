function [dx, dy, dz] = VectorDistanceROI(Node, ROI) % Вектор расстояния между двумя узлами
    dx = (Node(1)-ROI(:,1));
    dy = (Node(2)-ROI(:,2));
    dz = (Node(3)-ROI(:,3));
end