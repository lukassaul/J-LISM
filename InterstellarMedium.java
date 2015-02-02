/**
* This is an abstract class for modelling the interstellar medium
*
*  use this assuming that the distributions are homogenous
*  but not necessarily isotropic
*
*  extend this and change the distributions needed for the model
*/
public abstract class InterstellarMedium {

	public static double KB = 1.38*Math.pow(10,-23);
	public static double MP = 1.67*Math.pow(10,-27);
	public static double PI = Math.PI;

	/**
	* Make this maxwellian for now
	*
	*/

	public double dist(double v1, double v2, double v3) {
		return 0;
	}
	public double heliumDist(double v1, double v2, double v3) {
		return 0;
	}

	public double protonDist(double v1, double v2, double v3) {
		return 0;
	}

	public double hydrogenDist(double v1, double v2, double v3) {
		return 0;
	}

	public double deuteriumDist(double v1, double v2, double v3) {
		return 0;
	}

	public double electronDist(double v1, double v2, double v3) {
		return 0;
	}

	public double alphaDist(double v1, double v2, double v3) {
		return 0;
	}

	// etcetera...+

}
