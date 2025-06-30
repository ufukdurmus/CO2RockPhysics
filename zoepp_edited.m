function[pp,ps,ss,sh,thetap,thetas]=zoepp_edited(alpha1, beta1, rho1, alpha2, beta2, rho2)
%	--->>>   Calculates Zoeppritz reflection coefficients
%		   at an interface as a function of ray parameter
%
%		Assumed: iso / iso
%
%		from Aki-Richards, 1980, Chapter 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	--->>>  Velocities (p=alpha, s=beta), density=rho
%		 1: above
%		 2: below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	--->>>  specify maximum ray parameter
%	      ------------------------------------
	pmax = 1  / beta1;   % reflection coefficients will be complex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	--->>>    linearly interpolate to get ray parameters
%		----------------------------------------------
	p = linspace( 0.0, pmax, 100)';

	np = length(p);

	pp = zeros(np,1);  ps = zeros(np,1); sp = zeros(np,1); ss = zeros(np,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	--->>>   Loop over ray parameter,  get P-SV reflection coefficients
%		-----------------------------------------------------------
%		see Aki-Richards,  section 5.2.4
%		------------------------------------------
	for k = 1:np

	   [ rc ] = subZoeppritz( p(k), alpha1, beta1, rho1, alpha2, beta2, rho2 ) ;

	  pp(k) = rc(1,1);
	  ss(k) = rc(2,2);
	  ps(k) = rc(2,1);
	  sp(k) = rc(1,2);
	end

		sh = zeros( np,1 );

	for k = 1:np
	  	j1 = asin( p(k) * beta1 ) ;
		j2 = asin( p(k) * beta2 ) ;

		imp1 = rho1 * beta1 * cos( j1 ) ;
		imp2 = rho2 * beta2 * cos( j2 ) ;

		sh(k) = ( imp1 - imp2 ) / ( imp1 + imp2 ) ;
	end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	--->>>   reflection coefficients may be complex depending upon
%		  the maximum ray parameter
%	--->>>   plot the real part: ie  real(pp),  real(ss) ... etc
%
%	--->>> to display as function of incidence angle (rather than ray parameter)
%		thetap = asin( p * alpha1 ) * 180 / pi ;
%		thetas = asin( p *  beta1 ) * 180 / pi ;
%
		thetap = asin( p * alpha1 ) * 180 / pi ;
		thetas = asin( p *  beta1 ) * 180 / pi ;

end