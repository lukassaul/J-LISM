/**
*  Use this to compute the original velocity vector at infinity,
*    given a velocity and position at a current time.
*
*   // note velocity vector is reversed, represents origin to point at t=-infinity
*
*   Use static heliovector methods to avoid the dreaded "new Object()"...
*
*
*    !!!  tested extensively May 2004 !!!
*    not extensively enough.  Starting again - sept 2004
*/
private HelioVector getHPInfinity(HelioVector hpr, HelioVector hpv) {
	//System.out.println("test!");
	counter++;
	r = hpr.getR();
	if (r==0.0) r=Double.MIN_VALUE;

	// - define orbital plane as hpn
	// orbital plane is all vectors r with: r dot hpn = 0
	// i.e. hpn = n^hat (normal vector to orbital plane)
	hpn.setCoords(hpr);
	HelioVector.crossProduct(hpn, hpv); // defines orbital plane
	HelioVector.normalize(hpn);

	// We are going to set up coordinates in the orbital plane..
	// unit vectors in this system are hpox and hpoy
	// choose hpox so the particle is on the +x axis
	// such that hpox*hpn=0 (dot product)    theta of the given position is zero
	//hpox.setCoords(HelioVector.CARTESIAN, hpn.getY(), -hpn.getX(), 0);
	hpox.setCoords(HelioVector.CARTESIAN, hpr.getX()/r, hpr.getY()/r, hpr.getZ()/r);

	// this is the y axis in orbital plane
	hpoy.setCoords(hpn);
	HelioVector.crossProduct(hpoy, hpox);  // hopy (Y) = hpn(Z) cross hpox (X)


	// ok, we have defined the orbital plane !
	// now lets put the given coordinates into this plane
	// they are:  r, ptheta, rdot, thetadot
	//

	// we want rDot.  This is the component of traj. away from sun
	//
	rdot = hpv.dotProduct(hpr)/r;

	goingIn = true;
	if (rdot>=0) goingIn = false;


	// what is thetaDot?
	//
	//  use v_ = rdot rhat + r thetadot thetahat...
	//
	//  but we need the thetaHat ..
	//    thetaHat is just (0,1) in our plane
	//    because we are sitting on the X axis
	//    thats how we chose hpox
	//
	thetaHat.setCoords(hpoy);
	thetadot = hpv.dotProduct(thetaHat)/r;


	// NOW WE CAN DO THE KEPLERIAN MATH
	// let's calculate the orbit now, in terms of eccentricity e
	//  ,semi-major axis a, and rdoti (total energy = 1/2(rdoti^2))
	//
	//  the orbit equation is a conic section
	//  r(theta) = a(1-e^2) / [ 1 + e cos (theta - theta0) ]
	//
	//  thanks Bernoulli..
	//
	//   note:  a(1-e^2) = L^2/GM
	//
	//  	 rMin = a(1-e)
	//
	//			for hyperbolics, e>1
	//
	L=r*r*thetadot; // ok - we know our angular momentum (per unit mass)

	E=hpv.dotProduct(hpv)/2.0 - gmodMs/r;
	if (E <= 0) {
		d("discarding elipticals...");
		return null;
	}

	// speed at infinity better be more negative than current speed!
	if (rdoti > hpv.getR()) o("rdoti < hpv.getR() !!  ");

	if (thetadot==0.0) {
		rMin = 0.0;
		//o("Trajectory entirely radial! ");
	}
	else {
		rMin= (Math.sqrt(gmodMs*gmodMs/E/E+2.0*L*L/E)-gmodMs/E)/2.0;
		//rMin = Math.sqrt(2.0/thetadot/thetadot*(E+gmodMs/r));
	}
//	d("rmin:" + rMin);
	// the speed at infinity is all in the radial direction, rdoti (negative)
	rdoti=-Math.sqrt(2.0*E);




	// rMin had better be smaller than R..
	if (rMin > r) {
	//	d("rMin > r !! ");
		rMin = r;
	}

	// e and a now are available..
	e = L*L/rMin/gmodMs - 1.0;
	oneOverE = 1/e;
	if (e<=1) {
		d("didn't we throw out ellipticals already?? e= " + e);
		return null;
	}

	//a = rMin/(1.0-e);

	// do some debugging..
	//
	/*
	d("rdoti:  " + rdoti);d("r: " + r + "   L: " + L);
	d("ke: " + hpv.dotProduct(hpv)/2.0);
	d("pe: " + gmodMs/r);
	d("ke - pe: " + (hpv.dotProduct(hpv)/2.0 - gmodMs/r));
	d("rke: " + rdot*rdot/2.0);
	d("energy in L comp. " + L*L/r/r/2.0);
	d("ke - (rke+thetake): " + (hpv.dotProduct(hpv)/2.0 - rdot*rdot/2.0 - L*L/r/r/2.0));
	d("E dif: " + ( (rdot*rdot/2.0 + L*L/r/r/2.0 - gmodMs/r) - E ));
	d("E dif: " + ( (rdot*rdot/2.0 + thetadot*thetadot*r*r/2.0 - gmodMs/r) - E ));
	d("thetadot: " + thetadot);
	d("e: " + e);
	d("rMin: " + rMin+ "a: " + a + " e: " + e);
	*/


	// WE HAVE EVERYTHING in the orbit equation now except theta0
	//  the angle at which the orbit reaches rMin
	// now we use our current position to solve for theta - theta0
	//arg = a*(1.0-e*e)/e/r - 1.0/e;
	arg = rMin*(1.0+e)/e/r - oneOverE;
	// could also be
	arg = L*L/r/e/gmodMs - oneOverE;

	//d("arg: " + arg);

	theta0 = Math.acos(arg);
	// correct for going in when sitting really on x axis
	// we have the angle but we need to make sure we get the sign right
	// what with implementing the mulitvalued acos and all..
	//  acos returns a value from 0 to pi
	if (!goingIn) theta0 = -theta0;

	//d("theta0: "+ theta0);


	// the angle "swept out by the particle" to its current location is
	//
	//  this is also useful for the ionization calculation
	//
	//if (thetadot>0)  pthetai = theta0 + Math.acos(-oneOverE);
	//else  pthetai = theta0 - Math.acos(-oneOverE);

	//if (!goingIn) {
	//	angleSwept = Math.abs(pthetai - theta0) + Math.abs(Math.PI-theta0);
	//}
	//else {
	//	angleSwept = Math.abs(pthetai - Math.PI);
	//}
	//if (angleSwept>=Math.PI) o("angleSwept too big : " + angleSwept);
	angleSwept = Math.acos(-oneOverE);
	pthetai = -angleSwept + theta0;
	//if (!goingIn) pthetai = angleSwept - theta0;
//	d("angle swept: "+ angleSwept);
//	d("pthetai: " + pthetai);



	// now we know everything!!  The vector to the original v (t=-infinity)
	//
	// Vinf_ = - rdoti * rHat ( thetaInf )
	//
	vinfxo = rdoti*Math.cos(pthetai);
	vinfyo = rdoti*Math.sin(pthetai);

	// put together our answer from the coords in the orbital plane
	//  and the orbital plane unit vectors
	answer.setCoords(hpox);
	HelioVector.product(answer,vinfxo);
	HelioVector.product(hpoy,vinfyo);  // destroy hpoy here
	HelioVector.sum(answer,hpoy);
	//HelioVector.invert(answer);
	return answer;
}