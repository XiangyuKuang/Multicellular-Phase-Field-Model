function str = dot2d(num)
if ischar(num)
  str=num;
else
  str = num2str(num);
end
str(str == '.') = 'd';