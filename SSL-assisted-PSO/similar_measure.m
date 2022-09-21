function similar = similar_measure(x,y,d)
  sum_two = 0;
  sum_x = 0;
  sum_y = 0;
  sum_xx = 0;
  sum_yy = 0;
  for i=1:d
      sum_two = sum_two + x(i)*y(i);
      sum_x = sum_x + x(i);
      sum_y = sum_y + y(i);
      sum_xx = sum_xx + x(i)*x(i);
      sum_yy = sum_yy + y(i)*y(i);
  end
  similar = (sum_two - (sum_x * sum_y)/d)/(((sum_xx - sum_x^2/d)*(sum_yy - sum_y^2/d))^0.5);
end