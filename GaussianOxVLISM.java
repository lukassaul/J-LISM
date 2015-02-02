/**
*  The usual "hot model" here..  only with one species for now
*
*/
public class GaussianOxVLISM extends InterstellarMedium {

	public HelioVector vBulk1;
	public HelioVector vBulk2;
	public double temp1, temp2, density, fraction;
	public static double KB = 1.3807*Math.pow(10,-23);   // J/K
	public static double MP = 1.672621*Math.pow(10,-27);   // kg
	private double Mover2KbT1, Mover2KbT2, norm1, norm2, tt11, tt21, tt31, tt12, tt22, tt32, tt1, tt2;
	public static double NaN = Double.NaN;


	/**
	* Construct a gaussian VLISM
	*
	*  Here we put a primary and secondary component
	*/
	public GaussianOxVLISM(HelioVector _vBulk1, HelioVector _vBulk2, double _density, double _temp1, double _temp2, double _fraction) {

		vBulk1 = _vBulk1;
		temp1 = _temp1;
		vBulk2 = _vBulk2;
		temp2 = _temp2;
		density = _density;
		fraction = _fraction;
		//System.out.println("created gaussian vlism with ro: " + density + " temp: " + temp + " bulk: " + vBulk.getR());
		// normalization for maxwellian distribution
		norm1 = (1.0-fraction)*density*Math.pow(16.0*MP/2.0/Math.PI/KB/temp1,3.0/2.0);
		norm2 = fraction*density*Math.pow(16.0*MP/2.0/Math.PI/KB/temp2,3.0/2.0);
		// argument of exponential requires this factor
		Mover2KbT1 = 16.0*MP/2.0/KB/temp1;
		Mover2KbT2 = 16.0*MP/2.0/KB/temp2;

	}

	/**
	*
	* Pass in a heliovector if you feel like it.. slower?
	*   probably not..
	*/
	public double heliumDist(HelioVector v) {
		return dist(v.getX(), v.getY(), v.getZ());
	}

	public double heliumDist(double v1, double v2, double v3) {
		return dist(v1,v2,v3);
	}

	/**
	* default - pass in 3 doubles
	*/
	public double dist(double v1, double v2, double v3) {
		if (v1==NaN | v2==NaN | v3==NaN) return 0.0;
		tt11 = v1-vBulk1.getX(); tt21 = v2-vBulk1.getY(); tt31 = v3-vBulk1.getZ();
		tt12 = v1-vBulk2.getX(); tt22 = v2-vBulk2.getY(); tt32 = v3-vBulk2.getZ();
		tt1 = tt11*tt11+tt21*tt21+tt31*tt31;
		tt2 = tt12*tt12+tt22*tt22+tt32*tt32;
		return norm1*Math.exp(-Mover2KbT1*tt1) + norm2*Math.exp(-Mover2KbT2*tt2);
	}

	public double dist(HelioVector hv) {
		return dist(hv.getX(), hv.getY(), hv.getZ());
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

	/**
	*    Finally well tested, Aug. 2004
	*/
	public static final void main(String[] args) {

		/*final GaussianOxVLISM gvli = new GaussianOxVLISM(new HelioVector(HelioVector.CARTESIAN,
				-20000.0,0.0,0.0), 100.0, 6000.0);
		System.out.println("Created GVLISM" );

		//System.out.println("Maxwellian VLISM created, norm = " + norm);
		MultipleIntegration mi = new MultipleIntegration();
		//mi.setEpsilon(0.001);
		FunctionIII dist = new FunctionIII () {
			public double function3(double x, double y, double z) {
				return gvli.heliumDist(x, y, z);
			}
		};

		System.out.println("Test dist: " + dist.function3(20000,0,0));

		double[] limits = new double[6];
		limits[0]=mi.MINUS_INFINITY; limits[1]=mi.PLUS_INFINITY;
		limits[2]=mi.MINUS_INFINITY; limits[3]=mi.PLUS_INFINITY;
		limits[4]=mi.MINUS_INFINITY; limits[5]=mi.PLUS_INFINITY;
		double[] limits2 = new double[6];
		limits2[0]=-100000.0; limits2[1]=100000.0;
		limits2[2]=-100000.0; limits2[3]=100000.0;
		limits2[4]=-100000.0; limits2[5]=100000.0;
		double z = mi.integrate(dist,limits,128  );
		double z2 = mi.integrate(dist,limits2,128);

		System.out.println("Integral should equal density= " + z + "  " + z2);


		FunctionIII distZw = new FunctionIII () {
			public double function3(double r, double p, double t) {
				HelioVector hv = new HelioVector(HelioVector.SPHERICAL,r,p,t);
				return gvli.heliumDist(hv)*r*r*Math.sin(t);
			}
		};
		double[] limitsZw = {0.0,100000.0,0.0,2*Math.PI,0,Math.PI};
		double ans = mi.mcIntegrate(distZw, limitsZw, 200000);

		System.out.println("Integral should equal density= " + ans + "  old: " + z + "  " + z2);*/

	}

}
