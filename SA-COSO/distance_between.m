function dist = distance_between(x,y,d)
  distance = 0;
  for i=1:d
      distance = distance + (x(i) - y(i))^2;
  end
  dist = (distance)^0.5;
end