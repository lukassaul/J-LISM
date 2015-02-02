
import java.io.File;
import java.lang.Math;
import java.util.*;


/**
* This class computes the ionization percentage
* simple model with r^-2 dependence
*
*
*/
public class SurvivalProbability {

	public static double Ms = 2 * Math.pow(10,30);
	public static double G = 6.67 * Math.pow(10,-11);
	public static double AU = 1.5* Math.pow(10,11);
	public double v0, q, beta0, r0, r1, theta0, theta1, tempTheta1, thetaPrime, mu;
	// psi = angle swept out , p = Angular Momentum / mass, v0 = vinfinity of Moebius/Rucinski
	// see Wu & Judge 1979 APJ

	public SurvivalProbability(double lossRate1AU, double _mu) {
		beta0 = lossRate1AU;  //average loss rate at 1AU
	 	mu = _mu;
	}


	/**
	*  THis routine computes the fraction of ions that will be left given point / trajectory
	*  NOTE:  parameters of distribution must be given in coordinates
	*    relative to inflow direction
	*
	*   Taken from Judge&Wu 1979
	*/
	public double compute(double r, double theta, double vr, double vt) {
		//System.out.println("trying lossRate.compute");
		// Wu & Judge 1979 APJ -   - -

		// this stuff is to figure out theta prime, the angle swept out by
		// particle on way in from infinity

		q = Math.sqrt((1-mu)*G*Ms/r);
		v0 = Math.sqrt((vr*vr)+(vt*vt)-(2*q*q));
		r0 = r*vt*vt/(q*q);
		r1 = q*q*r/(v0*v0);

		theta0 = Math.PI/2 + Math.atan(-v0*vt/(q*q));  //  pi/2 < theta0 < pi
		if (((theta0)>Math.PI)|((theta0)<Math.PI/2)) {
			System.out.println("theta0 out of range!! (sp) " + theta0);
		}

		tempTheta1 = Math.atan(Math.sqrt(r0/r1 + 2*r0/r - r0*r0/(r*r)) / (r0/r - 1));
		//  but we want 0 < theta1 < pi ...
		if (tempTheta1 >= 0) theta1 = tempTheta1;
		else theta1 = tempTheta1 + Math.PI;
		if ((theta1>Math.PI)|(theta1<0)) System.out.println("theta1 out of range!! (sp)");


		if (vr > 0) {
			thetaPrime = theta0 + theta1;
		}

		else if (vr < 0) {
			thetaPrime = theta0 - theta1;
		}

		else if (vr == 0) {
 			thetaPrime = theta0;
		}

		if (thetaPrime<0) System.out.println("thetaPrime out of range!!");
	    //System.out.println("almost done sp.comute");

		// we have thetaPrime- now it's just an exponential
		return Math.exp(-beta0*AU*AU*thetaPrime/Math.abs(r*vt));
	}


	/**
	* This one gives the survival probablility, given the angle swept out in it's orbit
	*
	*   calculation of this angle must be done elsewhere in this routine.
	*
	*/
	public double compute(double thetaPrime, double L) {
		return 0;
		//return Math.exp(-bet0*AU*AU*tetaPrime/
	}



	/*
	* For testing only
	*/
	public static void main(String [] args) {
		SurvivalProbability sp = new SurvivalProbability(.1,0.0);
		System.out.println(""+sp.compute(2*AU, 2, 100000, 100000));
	}
}