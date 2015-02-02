import java.util.*;
//import cern.jet.math.Bessel;
//import com.imsl.math.Complex;

/**
*  This calculates the velocity distribution of neutral atoms in the heliosphere.
*  (according to the hot model or interstellar neutral dist. of choice)
*
*  Also computes losses to ionization - based on 1/r^2 rate
*
*  Uses an input distribution class to determine LISM f.
*
*  Lukas Saul - Warsaw 2000
*
*  Winter 2001 or so - got rid of Wu & Judge code
*  (still based on their calculation of course)
*
*  Last Modified - October 2002 - moved lism distribution to separate class
*   reviving - April 2004
*  Last used:
*  Sept 15. 2008
*/
public class NeutralDistribution {

	public static double KB = 1.38065 * Math.pow(10,-23);
	public static double Ms = 1.98892 * Math.pow(10,30);  // kg
	public static double G = 6.673 * Math.pow(10,-11);  // m^3/s^2/kg
	public static double AU = 1.49598* Math.pow(10,11); //meters
	public static double NaN = Double.NaN;
	public static double MAX = Double.MAX_VALUE;

	private double mu, beta, Gmod;
	public static int monte_num = 1000;


	// this gives our original neutral distribution in LISM
	private InterstellarMedium im;

	// these variables are for getHPInfinity
	private double rdot, rdoti, thetadot, L, ptheta, pthetai, thc, vdif;
	private double e, a, phi, shift;
	private HelioVector hpr, hpv, hpn, hpox, hpoy, hpvi, hpdif, hpThetaHat, answer;
	public boolean goingIn;

	// other objects and primitives
	private InterstellarMedium vlism;
	private double E, rMin, theta0;
	private HelioVector tempHV, thetaHat, temp2;
	private double testX, testY, arg, vinfxo, vinfyo, angleSwept, r;

	private double survivalFactor;
	public double gmodMs;
	public int counter;  // how many getHPInfinitys have been called?

	public boolean debug=true;
	private double oneOverE;

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
		Gmod = G*(1.0-mu);

		oneOverE = 0.0;

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
		temp2 = new HelioVector();
		goingIn=false;

		ptheta = Math.PI;

		gmodMs = Gmod*Ms;
		//if (gmodMs==0.0) gmodMs = Double.MIN_VALUE;

		counter = 0;
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

		// first make some reference vectors to help simplify geometryschool hou
		hpr.setCoords(HelioVector.SPHERICAL, r, theta, phi);
		// take the reverse of the velocity to help keep it straight
		hpv.setCoords(HelioVector.CARTESIAN, vx, vy, vz);
		return dist (hpr,hpv);
	}

	public double dist(HelioVector hpr_, HelioVector hpv_) {
		//System.out.println("test");
		hpr.setCoords(hpr_); hpv.setCoords(hpv_);
		double r = hpr.getR();

		hpvi = getHPInfinity(hpr, hpv); // here's the hard part - this call also sets thetaDif
		                               // the angle swept out by the particle from -inf to now

		if (hpvi==null) return 0.0;
		if ( !(hpvi.getX()<MAX && hpvi.getX()>-MAX) ) return 0.0;
		if ( !(hpvi.getY()<MAX && hpvi.getY()>-MAX) ) return 0.0;
		if ( !(hpvi.getZ()<MAX && hpvi.getZ()>-MAX) ) return 0.0;

		//System.out.println("made it here");
		// we need to multiply by the survival probability here...
		//System.out.println("thetadot: " + thetadot);
		//System.out.println("r: " + r);
		survivalFactor = Math.exp(-beta*AU*AU*angleSwept/L);
		//o("survival Factor: " + survivalFactor + " angleSwept: " + angleSwept + " L: " + L);

		// then do the gaussian with survival probability
		//double answer = vlism.heliumDist(hpvi.getX(),hpvi.getY(),hpvi.getZ())*survivalFactor;
		//o(""+answer);
		//System.out.println(hpvi.getX() + " " + hpvi.getY() + " " + hpvi.getZ());

		try {
			return vlism.dist(hpvi.getX(),hpvi.getY(),hpvi.getZ())*survivalFactor;
		}
		catch(Exception e) {
			return 0.0;
		}
	}



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


	/**
	* Use this to input a test for getHPInfinity..
	*
	*  That routine is driving me nuts.
	*
	*  Tested!  Working!  01:21 Dec4 2007
	*/
	public void testI(double x, double y, double z, double vx, double vy, double vz) {
		//HelioVector tttt = new HelioVector();
		HelioVector r1 = new HelioVector(HelioVector.CARTESIAN, x,y,z);
		HelioVector v1 = new HelioVector(HelioVector.CARTESIAN, vx,vy,vz);
		HelioVector tttt = getHPInfinity(r1,v1);
		o("r: " + r1.o());
		o("v: " + v1.o());
		o("goingIN: "+ goingIn + " theta0: " +theta0 + " pthetai: "+pthetai);
		o("angle swept: " + Math.acos(-1.0/e) + " rdoti: " + rdoti);
		if (tttt!=null) {
			o("v_inf : " +tttt.o());
			double einit= 0.5*v1.dotProduct(v1) - Ms*Gmod/Math.sqrt(r1.dotProduct(r1));
			double efinal= 0.5*tttt.dotProduct(tttt);
			o("energy difference = " + Math.abs(efinal-einit));
			o("\n\n\n");
		}
	}


	/**
	*  Here is the formula to graph
	*  r(theta) = a(1-e^2) / [ 1 + e cos (theta - theta0)
	*
	*/
	public void graphOrbit(HelioVector v_inf) {	// not implemented
	}


	/**
	*  Use this to do a more complete test of getHPInfinity
	*/
	/*public void testIII() {
		// first let's try a sweep in position phi coord, leaving v constant.
		//  do a hundred for fun
		file fu = new file("testIII.dat");
		fu.initWrite(false);
		for (int i=0; i<100; i++) {
			HelioVector r1 = new HelioVector(HelioVector.SPHERICAL, AU,i*2*Math.PI/100,0.0);
			HelioVector v1 = new HelioVector(HelioVector.CARTESIAN, 30000.0,0.0,0.0);
			HelioVector tttt = getHPInfinity(r1,v1);
			fu.write(i*2*Math.PI/100 + "\t" +
	*/


	/**
	* for testing only
	*/
	public static final void main(String[] args) {

		//	o("MAX_VALUE: " + Double.MAX_VALUE);
		double test234=3.2/0.0;
		Date dt3 = new Date();
		MultipleIntegration mi = new MultipleIntegration();


		// **********************************
		// generate 2D color countour plot of heliospheric density of neutrals
		//		//
		file f = new file ("norm_test1_6500_26000.txt");
		f.initWrite(false);
		int res = 20;
		for (double xx=-5.0*AU; xx<=5.0*AU; xx+= 10.0*AU/res) {
			for (double yy=-5.0*AU; yy<=5.0*AU; yy+= 10.0*AU/res) {
				final double x = xx;
				final double y = yy;

				FunctionIII f3 = new FunctionIII() {
					public double function3(double a, double b, double c) {

						HelioVector bulk = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
						NeutralDistribution ndH = new NeutralDistribution(new GaussianVLISM(bulk,0.015,6500.0),  0.0, 1.0*Math.pow(10,-7));
						ndH.debug=false;
						double tbr = ndH.dist(new HelioVector(HelioVector.CARTESIAN, x,y,0.0),
										new HelioVector(HelioVector.CARTESIAN,a,b,c));
						//System.out.println(tbr);
						return tbr;
					}
				};

				//System.out.println("test: " + f3.function3(500000,500000,500000));
				//HelioVector testPos;
				//HelioVector testSpeed;
				//double minf = MultipleIntegration.MINUS_INFINITY;
				//double pinf = MultipleIntegration.PLUS_INFINITY;

				//MultipleIntegration mi = new MultipleIntegration();

				//double[] limits = {minf,pinf,minf,pinf,minf,pinf};
				//for (double lim=50000.0; lim<110000.0; lim+=10000.0) {
				//	double[] limits2 = {-lim,lim,-lim,lim,-lim,lim};
				//	System.out.println("lim: " + lim + " ans: " +mi.integrate(f3,limits2,32));
				//}
				// let's go with 70000 for now at 32 steps

				double lim = 80000.0;
				double[] limits2 = {-lim,lim,-lim,lim,-lim,lim};

				//System.out.println("lim: " + lim + " ans: " +mi.integrate(f3,limits2,32));

				double zz =mi.mcIntegrate(f3,limits2,monte_num);
				System.out.println("done x: "+xx+" y: "+yy+" dens: "+zz);
				f.write(xx+"\t"+yy+"\t"+zz+"\n");
			}
		}
		f.closeWrite();
		// *******************************


		/*
		HelioVector bulk = new HelioVector(HelioVector.CARTESIAN,26000.0, 0.0, 0.0);
		final HelioVector posVec = new HelioVector(HelioVector.SPHERICAL, 300*AU, 0.0, 0.0);
		final NeutralDistribution nd = new NeutralDistribution(new GaussianVLISM(bulk,0.015,6300.0),  0.8, 1.1*Math.pow(10,-7));
		nd.debug=false;
		*/


		//	o("2 / min : " + 2.0/Double.MIN_VALUE);
		// bulk speed
		//HelioVector bulk = new HelioVector(HelioVector.SPHERICAL,26000.0,74.5*Math.PI/180.0,-5.0*Math.PI/180);
		// let's make three of these, H He O
		// set up the interstellar distribution
		//	NeutralDistribution ndH = new NeutralDistribution(new GaussianVLISM(bulk,0.19,6300.0),  1.0, 7.4*Math.pow(10,-7));

		//final NeutralDistribution nd = new NeutralDistribution(new KappaVLISM(bulk,0.015,6300.0, 1.6),  0.0, 1.1*Math.pow(10,-7));

		// going to look at particles at a few energies in Earth frame and see their total deflection
		/*
		double MP = 1.672621E-27;
		double EV = 1.60217646E-19;


		double ev1000v = Math.sqrt(2*1000/MP*EV);
		System.out.println("should be 440km/s: " + ev1000v);

		double[] energies = {14,20,30,40,50,75,100,125,150,200};
		double[] speeds = new double[energies.length];

		for (int i=0; i<speeds.length; i++) {
			speeds[i] = Math.sqrt(2*energies[i]/MP*EV) - 30000.0;
			System.out.println(speeds[i]);
			HelioVector velVecta = new HelioVector(HelioVector.CARTESIAN, 0.0, speeds[i], 0.0);
			HelioVector posVecta = new HelioVector(HelioVector.CARTESIAN, 1.0*AU, 0.0, 0.0);

			HelioVector infVecta = nd.getHPInfinity(posVecta, velVecta);
			//System.out.println(infVecta);

			if (infVecta!= null) System.out.println(energies[i]+ "\t" + (180.0/Math.PI*velVecta.getAngle(infVecta)));
		}
		*/



		//NeutralDistribution ndO = new NeutralDistribution(new GaussianVLISM(bulk,0.00005,6300.0),  0.0, 7.5*Math.pow(10,-7));
		//ndH.debug=false;
		//ndHe.debug=false;
		//ndO.debug=false;


		//TEST SPHERICAL LIMITS
		/*double[] limits = { 0.0, 300000.0, 0.0, 2.0*Math.PI, 0.0, Math.PI };
		//o("doing pos: " + pos);
		//final HelioVector posVec = new HelioVector(HelioVector.SPHERICAL, 200*AU, 0.0, 0.0);
		FunctionIII f3 = new FunctionIII() {
			public double function3(double v, double p, double t) {
				double tbr = ndHe.dist(posVec, new HelioVector(HelioVector.SPHERICAL,v,p,t))*v*v*Math.sin(t);
				// that is the integrand of our v-space integration
				return tbr;
			}
		};
		o( mi.mcIntegrate(f3, limits, 10000) + "\n");
		*/




		// Determine distribution of speeds at a point
		/*file newF = new file("speed_distribution_2.txt");
		newF.initWrite(false);
		HelioVector bulk = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		final NeutralDistribution ndH = new NeutralDistribution(new GaussianVLISM(bulk,0.015,6300.0),  0.0, 1.0*Math.pow(10,-7));
		ndH.debug=false;
		//
		// set position
		final double x_pos = 100.0*AU;
		final double y_pos = 100.0*AU;
		final double z_pos = 100.0*AU;
		MultipleIntegration mi = new MultipleIntegration();double lim = 80000.0;
		double[] limits2 = {-lim,lim,-lim,lim};
		for (double v=0; v<80000; v+=(80000.0/20.0)) {  // step through speeds
			// at each speed go through
			final double ttv = v;
			FunctionII f2 = new FunctionII() {
				public double function2(double b, double c) {
					double tbr = ndH.dist(new HelioVector(HelioVector.CARTESIAN, x_pos, y_pos, z_pos),
							new HelioVector(HelioVector.CARTESIAN,Math.sqrt(ttv*ttv-b*b-c*c),b,c));
					//System.out.println(tbr);
					return tbr;
				}
			};
			newF.write(v+"\t"+mi.mcIntegrate(f2,limits2,10000)+"\n");
		}
		newF.closeWrite();*/


	    // COMPARE FLUX AT POINTS AROUND SUN
		/*file f = new file("testrestults.txt");
		f.initWrite(false);
		for (double phi=0; phi<2.0*Math.PI; phi+=8*Math.PI/360.0) {
			//o("trying phi: " + phi);
			testPos = new HelioVector(HelioVector.CARTESIAN, AU*Math.cos(200*Math.PI/180), AU*Math.sin(200*Math.PI/180),0.0);
			testSpeed = new HelioVector(HelioVector.CARTESIAN, 50000*Math.cos(200*3.14/180+Math.PI/2.0)*Math.cos(phi), 50000*Math.sin(200*3.14/180+Math.PI/2)*Math.cos(phi), 50000*Math.sin(phi));
			f.write(phi*180/3.14159+"\t"+ndH.dist(testPos,testSpeed)+"\t"+ndHe.dist(testPos,testSpeed)+"\n");
		}
		f.closeWrite();
		*/


		/*
		file f = new file("testrestults.txt");
		f.initWrite(false);
		for (double phi=0; phi<2.0*Math.PI; phi+=Math.PI/200.0) {
			//o("trying phi: " + phi);
			testPos = new HelioVector(HelioVector.CARTESIAN, AU,0.0,0.0);
			testSpeed = new HelioVector(HelioVector.CARTESIAN, 49000*Math.cos(phi+Math.PI/2.0), 49000*Math.sin(phi+Math.PI/2), 0.0);
			HelioVector hpinf = ndHe.getHPInfinity(testPos,testSpeed);
			if (hpinf!=null) f.write(phi*180/3.14159+"\t"+hpinf.getR()+"\t"+bulk.dotProduct(hpinf)+"\n");
		}
		f.closeWrite();
		*/


  		//  TEST OF GET HP INFINITY BELOW
	/*	Date dt4 = new Date();
		o("constructor took: "+(dt4.getTime()-dt3.getTime()));

		System.out.println(nd.dist(0.0,0.0,0.0,0.0,0.0,0.0));

		Date dt1 = new Date();

		nd.testI(0.0,1*AU,0.0,  -10000,0.0,0.0);
		nd.testI(0.0,1.0,0.0,  10000,0.0,0.0);
		nd.testI(0.0,1.0,0.0,  0.0,10000,0.0);
		nd.testI(0.0,1.0,0.0,  0.0,-10000,0.0);
		nd.testI(0.0,1.0,0.0,  0.0,0.0,10000);
		nd.testI(0.0,1.0,0.0,  0.0,0.0,-10000);
	*/


		// let's test going straight in from six cardinal directions..
/*
		nd.testI(0.0,AU,0.0,  0.0,-20000.0,200000.0);
		nd.testI(0.0,AU,0.0,  0.0,20000.0,200000.0);
		nd.testI(0.0,-AU,0.0, 0.0,-20000.0,200000.0);
		nd.testI(0.0,-AU,0.0, 0.0,20000.0,200000.0);
		nd.testI(-AU,0.0,0.0,  0,0.-200000.0,0.0);
		nd.testI(-AU,0.0,0.0,  0.0,0.0,200000.0);
		nd.testI(-AU,0.0,0.0,  0.0,0.0,-200000.0);
		nd.testI(0.0,-AU,0.0,  0.0,50000.0,0.0);
		nd.testI(0.0,AU,0.0,   0.0,-50000.0,0.0);
		nd.testI(0.0,0.0,-AU,  0.0,0.0,50000.0);
		nd.testI(0.0,0.0,AU,   0.0,0.0,-50000.0);


		nd.testI(0.0,-AU,0.0,  0.0,0.0,50000);
		nd.testI(0.0,-AU,0.0,  0.0,0.0,100000);
		nd.testI(0.0,-AU,0.0,  0.0,0.0,500000);

		nd.testI(-AU,0.0,0.0,  0.0,0.0,50000);
		nd.testI(-AU,0.0,0.0,  0.0,0.0,100000);
		nd.testI(-AU,0.0,0.0,  0.0,0.0,500000);
*/

		//	Date dt2 = new Date();
		//	System.out.println("6 took: "+(dt2.getTime()-dt1.getTime()));




		// we are going to generate some line plots..
		/*
		int ts = 10;
		double[] y = new double[ts];
		double[] x = new double[ts];

		// a point at 1AU for starters
		HelioVector pos = new HelioVector(HelioVector.CARTESIAN,AU,0,0);
		for (int i=0; i<y.length; i++) {
			HelioVector vel = new HelioVector(HelioVector.CARTESIAN,i*1000,0,0);
			x[i]=i*100;
			y[i]=nd.dist(pos,vel);
			System.out.println(x[i]+"\t"+y[i]);
		}*/

	}

	private static void o(String s) {
		System.out.println(s);
	}

	private void d(String s) {
			if (debug) System.out.println(s);
	}
}