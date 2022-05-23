function [number, idx, distances, minDist] = findClosestValue(value, array)

	distances    = abs(array - value);
	minDist = min(distances);
	idx     = find(distances == minDist);

	number = array(idx);
	

end