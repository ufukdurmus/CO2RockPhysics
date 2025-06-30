	function [ rc ] = subZoeppritz( p, a1, b1, r1, a2, b2, r2 )

		i1 = asin( p * a1 ) ;
		j1 = asin( p * b1 ) ;
	 	i2 = asin( p * a2 ) ;
		j2 = asin( p * b2 ) ;

	M(1,:) = [  -a1*p   -cos(j1)    a2*p    cos(j2) ];
	M(2,:) = [ cos(i1)  -b1*p     cos(i2)    -b2*p   ];

	M(3,1) = 2*r1*b1*b1*p*cos(i1);
	M(3,2) =  r1*b1*(1-2*b1*b1*p*p) ;
	M(3,3) =  2*r2*b2*b2*p*cos(i2)  ;
	M(3,4) = r2*b2*(1-2*b2*b2*p*p)  ;
	
	M(4,1) = -r1*a1*(1-2*b1*b1*p*p) ;
	M(4,2) = 2*r1*b1*b1*p*cos(j1) ;
	M(4,3) = r2*a2*(1-2*b2*b2*p*p) ;
	M(4,4) = -2*r2*b2*b2*p*cos(j2) ;

	N(1,:) = [  a1*p    cos(j1)    -a2*p    -cos(j2) ];
	N(2,:) = [ cos(i1)   -b1*p     cos(i2)  -b2*p    ];

	N(3,1) = 2*r1*b1*b1*p*cos(i1) ;
	N(3,2) = r1*b1*(1-2*b1*b1*p*p);
	N(3,3) = 2*r2*b2*b2*p*cos(i2);
	N(3,4) = r2*b2*(1-2*b2*b2*p*p) ;
	
	N(4,1) = r1*a1*(1-2*b1*b1*p*p) ;
	N(4,2) = -2*r1*b1*b1*p*cos(j1);
	N(4,3) = -r2*a2*(1-2*b2*b2*p*p);
	N(4,4) = 2*r2*b2*b2*p*cos(j2) ;

	rc = M \ N ;

	return
	 
