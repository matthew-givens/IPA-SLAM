function vX = X( v )
%  Calculates the skew-symmetrix cross product matrix from a vector

vX = [   0 , -v(3),  v(2);
	   v(3)     0 , -v(1);
	  -v(2)   v(1),    0];

end

