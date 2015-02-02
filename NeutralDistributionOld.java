import java.util.*;
//import cern.jet.math.Bessel;
//import com.imsl.math.Complex;

/**
*  This calculates the velocity distribution of neutral atoms in the heliosphere.
*  (according to the hot model or interstellar neutral dist. of choice)
*
*  Also computes losses to ionization!
*
*  Uses an input distribution class to determine LISM f.
*
*  Lukas Saul - Warsaw 2000
*
*  Winter 2001 or so - got rid of Wu & Judge code
*  Last Modified - October 2002 - moved lism distribution to separate class
*   reviving - April 2004
*/
public class NeutralDistribution {

	public static double Ms = 1.98892 * Math.pow(10,30);  // kg
	public static double G = 6.673 * Math.pow(10,-11);
	public static double AU = 1.49598* Math.pow(10,11); //meters
	public static double NaN = Double.NaN;

	private double mu, beta, Gmod;


	// this gives our original neutral distribution in LISM
	private InterstellarMedium im;

	// these variables are for getHPInfinity
	private double rdot, rdoti, thetadot, L, ptheta, pthetai, thc, vdif;
	private double e, a, phi, shift;
	private HelioVector hpr, hpv, hpn, hpox, hpoy, hpvi, hpdif, hpThetaHat, answer;
	public boolean goingIn;


	// other objects and primitives
	private InterstellarMedium vlism;
	private double test2;
	private HelioVector tempHV, thetaHat;
	private double testX, testY, arg;


	/**
	*  Parameters:
	*    InterstellarMedium (starting point)
	*    radiationToGravityRatio (mu)
	*    SurvivalProbability (ionization rate calculatore)
	*
	*/
	public NeutralDistribution(InterstellarMedium im, double radiationToGravityRatio, double beta0) {

		mu = radiationToGravityRatio;
		beta = beta0;

		//sp = new SurvivalProbability(beta, mu);
		vlism = im;

		Gmod = G*(1-mu);



		// initialize them now.  We don't want to build java objects during calculations.
		hpr = new HelioVector();
		hpv = new HelioVector();
		hpn = new HelioVector();
		hpox = new HelioVector();
		hpoy = new HelioVector();
		hpvi = new HelioVector();
		hpdif = new HelioVector();
		hpThetaHat = new HelioVector();
		answer = new HelioVector();
		tempHV = new HelioVector();
		thetaHat = new HelioVector();
		goingIn=false;
	}



	/**
	*  Calculate the full distribution.  This time perform all transformations
	*   explicitly, instead of reading them from a paper (WU and JUDGE).
	*   Coordinates are expected in Heliospheric coordinates (phi=0,theta=0 = direct upwind)
	*
	*   Uses HelioVector.java to assist in geometry calculations
	*
	*   we are being careful now not to create "new" variables (instances) to save time
	*     especially our vectors (HelioVectors)
	*/
	public double dist(double r, double theta, double phi, double vx, double vy, double vz) {
		//o("Calling dist: "+r+" "+theta+" "+phi+" "+vx+" "+vy+" "+vz);

		// first make some reference vectors to help simplify geometry
		hpr.setCoords(HelioVector.SPHERICAL, r, theta, phi);
		// take the reverse of the velocity to help keep it straight
		hpv.setCoords(HelioVector.CARTESIAN, vx, vy, vz);

		hpvi = getHPInfinity(hpr, hpv); // here's the hard part - this call also sets thetaDif
		                               // the angle swept out by the particle from -inf to now

		if (hpvi==null) return 0;

		// then calculate the v0 - vinf
		//hpdif.setCoords(hpInflow);
		//HelioVector.difference(hpdif, hpvi);
		//o("hpdif: "+hpdif);

		// we need to multiply by the survival probability here...

		// then do the gaussian with survival probability
		double answer = vlism.heliumDist(hpvi.getX(),hpvi.getY(),hpvi.getZ());
		//o(""+answer);
		return answer;
	}



	/**
	*  Use this to compute the original velocity vector at infinity,
 	*    given a velocity and position at a current time.
	*
	*   // note velocity vector is reversed, represents origin to point at t=-infinity
	*
	*   Use static heliovector methods to avoid the dreaded "new Object()"...
	*

	*/
	private HelioVector getHPInfinity(HelioVector hpr, HelioVector hpv) {
		double r = hpr.getR();

		// - define orbital plane
		// orbital plane is all vectors r with: r dot hpn = 0
		// i.e. hpn = n^hat (normal vector to orbital plane)
		hpn.setCoords(hpr);
		HelioVector.crossProduct(hpn, hpv); // defines orbital plane
		if (hpn.getR()==0) {return null;}  // problems with div. by zero if we are at the origin
		else HelioVector.normalize(hpn);


		// are we coming into the sun or going away?
		// we need to know this for the geometry;
		goingIn = false;
		if (hpr.dotProduct(hpv)<0) goingIn = true;


		// We are going to set up coordinates in the orbital plane..
		// unit vectors in this system are hpox and hpoy
		// choose hpox so the particle is on the -x axis
		// such that hpox*hpn=0 (dot product)
		//hpox.setCoords(HelioVector.CARTESIAN, hpn.getY(), -hpn.getX(), 0);
		hpox.setCoords(HelioVector.CARTESIAN, -hpr.getX()/r, -hpr.getY()/r, -hpr.getZ()/r);
		HelioVector.normalize(hpox);

		// this is the y axis in orbital plane
		hpoy.setCoords(hpn);
		HelioVector.crossProduct(hpoy, hpox);  // hopy (Y) = hpn(Z) cross hpox (X)


		o("hpoy size: " + hopy.getR());
		HelioVector.normalize(hpoy);
		o("hopy size: " + hopy.getR());

		// we want the ++ quadrant to be where the particle came from, i.e.dtheta/dt positive
		/*testY = hpr.dotProduct(hpoy);
		testX = hpr.dotProduct(hpox);
		if (goingIn && (testY<0)) {
			HelioVector.invert(hpoy);	//o("we have inverted hpoy 1...");
		}
		else if (!goingIn && (testY>0)) {
			HelioVector.invert(hpoy);	//o("we have inverted hpoy 2...");
		}
		if (goingIn && (testX<0)) {
			HelioVector.invert(hpox);	//o("we have inverted hpox 1...");
		}
 		else if (!goingIn &&(testY>0)) {
			HelioVector.invert(hpoy);	//o("we have inverted hpox 2...");
		}*/

		// test so far
		o("hpr:  " + hpr.o());
		o("hpv: " + hpv.o());
		o("hpox: " + hpox.o());
		o("hpoy: " + hpoy.o());

		o("hpr dot hpox: " + hpr.dotProduct(hpox));

		// ok, we have defined the orbital plane !
		// now lets put the given coordinates into this plane
		// they are:  r, ptheta, rdot, thetadot

		// Now we know ptheta, the orbital plane theta component
		//ptheta = hpr.getAngle(hpox);
		//if (!goingIn) ptheta = -ptheta + 2*Math.PI;

		ptheta = Math.PI;

		// we want rDot.  This is the component of traj. parallel to our position
		//
		rdot = -hpv.dotProduct(hpr)/r;


		// one way to get thetadot:
		//
		//   thetadot = v_ cross r_ / r^2    ????
		//
		// using a temp variable to avoid "new Object()..."
		//
		//tempHV.setCoords(hpv);
		//HelioVector.crossProduct(tempHV,hpr);
		//thetadot = tempHV.getR()/r/r;

		// yet a fucking nother way to get thetadot
		//
		//  is from v = rdot rhat + r thetadot thetahat...
		//
		//  thetadot thetahat = (v_ - rdot rhat) / r
		//
		//  but we need the thetaHat to get the sign right
		//    thetaHat is just (0,-1) in our plane
		//    because we are sitting on the minus X axis
		//    thats how we chose hpox
		//
		thetaHat.setCoords(hpoy);
		thetaHat.invert();   //we have thetahat

		tempHV.setCoords(hpv);
		HelioVector.product(tempHV,1.0/r);
		HelioVector.difference(tempHV,hpr.product(rdot/r));
		thetadot = tempHV.getR()/r;
		// but what is the sign?
		if (hpv.dotProduct(thetaHat),0) thetadot = -thetadot;

		// test so far
		o("Coords in orbital plane - r: " + hpr.getR());
		o("ptheta: " + ptheta);
		o("thetadot: " + thetadot);
		o("rdot: " + rdot);


		// NOW WE CAN DO THE KEPLERIAN MATH
		// let's calculate the orbit now, in terms of eccentricity e
		//  ,semi-major axis a, and rdoti (total energy = 1/2(rdoti^2))
		//
		//  the orbit equation is a conic section
		//  r(theta) = a(1-e^2) / [ 1 + e cos (theta - theta0) ]
		//
		//  thanks Bernoulli..

		L=r*r*thetadot; // ok - we know our angular momentum (/Mhe)


		//  NOW use the ENERGY conservation to compute
		// rdoti - the velocity at infinity
		//
		//if ((test2=rdot*rdot - L*L/r/r + 2*Gmod*Ms/r) < 0) {
		test2=hpv.dotProduct(hpv)-2*Gmod*Ms/r;
		o("E= " + test2);
		if (test2 < 0) {
			o("discarding elipticals...");
		//	System.out.println("discarding elipticals");
			return null;
		}


		o("E dif: " + ( (rdot*rdot + L*L/r/r - 2*Gmod*Ms/r) - test2 ));

		rdoti=-Math.sqrt(test2);

		o("rdoti: " + rdoti);


		a = -Gmod*Ms/rdoti/rdoti;
	//	e = Math.sqrt(1+L*L*rdoti*rdoti/Gmod/Gmod/Ms/Ms);
		e = Math.sqrt(1-L*L/a/Gmod/Ms);

		if (e<=1) {
			System.out.println("I thought we threw out ellipticals already??");
			return null;
		}

		double maxPhi = Math.acos(1/e);
		o("maxPhi: " + maxPhi);
		o("e= " + e + " a= " + a);
		// that pretty much defines the orbit, we need to find the angular shift
		// to match with our coordinate system

		arg = a*(1-e*e)/e/r - 1/e;
		o("arg: " + arg);
		if (arg<-1 | arg>1) {
			o("problems here 9a8nv3");
			//return null;
		}
		phi = Math.acos(arg)-maxPhi;

		if (goingIn) { // we are in the +y half, ptheta<180
			shift = ptheta-(Math.PI-phi);
			pthetai = (Math.PI-maxPhi)+shift;
		}
		else {
			shift = (Math.PI-phi)-(2*Math.PI-ptheta);
			pthetai = (Math.PI-maxPhi)+shift;
		}

		o("phi: " + phi + " pthetai" + pthetai);

		// now we know vinf in it's entirety.  The vector, in original coordinates:
		double vinfxo = rdoti*Math.cos(pthetai);
		double vinfyo = rdoti*Math.sin(pthetai);

		// put together our answer from the coords in the orbital plane
		//  and the orbital plane unit vectors
		answer.setCoords(hpox);
		HelioVector.product(answer,vinfxo);
		//HelioVector.sum(answer,hpoy.product(vinfyo));
		HelioVector.product(hpoy,vinfyo);
		HelioVector.sum(answer,hpoy);

		// that's a vector going away from the sun.  Invert it
		HelioVector.invert(answer);
		return answer;
	}

	/**
	* for testing only
	*/
	public static final void main(String[] args) {

		/*parameters:
		*/


		Date dt3 = new Date();

		NeutralDistribution nd = new NeutralDistribution(new GaussianVLISM(new HelioVector(),0.0,0.0),0.0,0.0);
		Date dt4 = new Date();
		o("constructor took: "+(dt4.getTime()-dt3.getTime()));

		// test getHPInfinity
		//HelioVector tttt = new HelioVector();
		Date dt1 = new Date();
		HelioVector r1 = new HelioVector(HelioVector.CARTESIAN, 0.01,AU, 0.01);
		HelioVector v1 = new HelioVector(HelioVector.CARTESIAN, 0.01, 0.01, -50000);
		o("Trying r = AUy, v = -50000x");
		HelioVector tttt = nd.getHPInfinity(r1,v1);
		o("goingIN: "+ nd.goingIn);
		o(tttt+"");
		double einit= 0.5*v1.dotProduct(v1) - Ms*nd.Gmod/Math.sqrt(r1.dotProduct(r1));
		double efinal= 0.5*tttt.dotProduct(tttt);
		o("energy difference = " + Math.abs(efinal-einit));
		o("\n\n\n");


		r1 = new HelioVector(HelioVector.CARTESIAN, AU,0.01,0.01);
		v1 = new HelioVector(HelioVector.CARTESIAN, 0.01,-50000,  0.01);
		o("Trying x axis, -50000y");
		tttt = nd.getHPInfinity(r1,v1);
		o("goingIN: "+ nd.goingIn);
		o(tttt+"");
		einit= 0.5*v1.dotProduct(v1) - Ms*nd.Gmod/Math.sqrt(r1.dotProduct(r1));
		efinal= 0.5*tttt.dotProduct(tttt);
		o("energy difference = " + Math.abs(efinal-einit));
		o("\n\n\n");


		r1 = new HelioVector(HelioVector.CARTESIAN, 0.01,0.01,AU);
		v1 = new HelioVector(HelioVector.CARTESIAN, 0.01, 0.01, -50000);
		o("Trying z axis, -50000z");
		tttt = nd.getHPInfinity(r1,v1);
		o("goingIN: "+ nd.goingIn);
		o(tttt+"");
		einit= 0.5*v1.dotProduct(v1) - Ms*nd.Gmod/Math.sqrt(r1.dotProduct(r1));
		efinal= 0.5*tttt.dotProduct(tttt);
		o("energy difference = " + Math.abs(efinal-einit));
		o("\n\n\n");


		r1 = new HelioVector(HelioVector.CARTESIAN, AU,0.01,0.01);
		v1 = new HelioVector(HelioVector.CARTESIAN, 50000, 0.01, 0.01);
		o("Trying x axis, 50000x");
		tttt = nd.getHPInfinity(r1,v1);
		o("goingIN: "+ nd.goingIn);
		o(tttt+"");
		einit= 0.5*v1.dotProduct(v1) - Ms*nd.Gmod/Math.sqrt(r1.dotProduct(r1));
		efinal= 0.5*tttt.dotProduct(tttt);
		o("energy difference = " + Math.abs(efinal-einit));
		o("\n\n\n");


		r1 = new HelioVector(HelioVector.CARTESIAN, 0.01,0.01,AU);
		v1 = new HelioVector(HelioVector.CARTESIAN, 0.01, 0.01, 50000);
		o("Trying z axis, 50000z");
		tttt = nd.getHPInfinity(r1,v1);
		o("goingIN: "+ nd.goingIn);
		o(tttt+"");
		einit= 0.5*v1.dotProduct(v1) - Ms*nd.Gmod/Math.sqrt(r1.dotProduct(r1));
		efinal= 0.5*tttt.dotProduct(tttt);
		o("energy difference = " + Math.abs(efinal-einit));
		o("\n\n\n");

		r1 = new HelioVector(HelioVector.CARTESIAN, 0.01,AU,0.01);
		v1 = new HelioVector(HelioVector.CARTESIAN, 0.01, 50000, 0.01);
		o("Trying y axis, 50000y");
		tttt = nd.getHPInfinity(r1,v1);
		o("goingIN: "+ nd.goingIn);
		o(tttt+"");
		einit= 0.5*v1.dotProduct(v1) - Ms*nd.Gmod/Math.sqrt(r1.dotProduct(r1));
		efinal= 0.5*tttt.dotProduct(tttt);
		o("energy difference = " + Math.abs(efinal-einit));
		o("\n\n\n");


		//OK, that's 6 easy tests to see whats happening.
		// to be really sure lets give her 3 more interesting cases and

		r1 = new HelioVector(HelioVector.CARTESIAN, 0.01,0.01,AU);
		v1 = new HelioVector(HelioVector.CARTESIAN, 0.01, 0.01, 0.01);
		o("Trying z axis, 0v");
		tttt = nd.getHPInfinity(r1,v1);
		o("goingIN: "+ nd.goingIn);
		o(tttt+"");
		einit= 0.5*v1.dotProduct(v1) - Ms*nd.Gmod/Math.sqrt(r1.dotProduct(r1));
		efinal= 0.5*tttt.dotProduct(tttt);
		o("energy difference = " + Math.abs(efinal-einit));
		o("\n\n\n");

		r1 = new HelioVector(HelioVector.CARTESIAN, 0.01,0.01,AU);
		v1 = new HelioVector(HelioVector.CARTESIAN, 50000, 0.01, 0.01);
		o("Trying z axis, 50000x");
		tttt = nd.getHPInfinity(r1,v1);
		o("goingIN: "+ nd.goingIn);
		o(tttt+"");
		einit= 0.5*v1.dotProduct(v1) - Ms*nd.Gmod/Math.sqrt(r1.dotProduct(r1));
		efinal= 0.5*tttt.dotProduct(tttt);
		o("energy difference = " + Math.abs(efinal-einit));
		o("\n\n\n");

		r1 = new HelioVector(HelioVector.CARTESIAN, 0.01,0.01,AU);
		v1 = new HelioVector(HelioVector.CARTESIAN, -50000, 0.01, 0.01);
		o("Trying z axis, -50000x");
		tttt = nd.getHPInfinity(r1,v1);
		o("goingIN: "+ nd.goingIn);
		o(tttt+"");
		einit= 0.5*v1.dotProduct(v1) - Ms*nd.Gmod/Math.sqrt(r1.dotProduct(r1));
		efinal= 0.5*tttt.dotProduct(tttt);
		o("energy difference = " + Math.abs(efinal-einit));
		o("\n\n\n");


		Date dt2 = new Date();
		System.out.println("9 took: "+(dt2.getTime()-dt1.getTime()));


// hello
	/*

		Date d0= new Date();
		o("r on z axis, v on -x");
		double d = nd.dist(1*AU,0,.001,0,0,0);

		o(d+"\n\n\n\n\nr on x axis, v on -x");
		double e = nd.dist(100*AU,0,Math.PI/2, 0,0,.001);

		double f = nd.dist(100*AU,0,.001, 0,0,.001);
		Date d1 = new Date();
		o(e+"\n"+f+"\n3 dists took: " + (d1.getTime()-d0.getTime()) );


		double dd = nd.density(AU,0,0);
		o(dd+"");
		Date d2 = new Date();
		o((d2.getTime()-d1.getTime()));


		*/

	}

	private static void o(String s) {
		System.out.println(s);
	}
}