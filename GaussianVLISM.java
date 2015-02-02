import java.util.Date;

/**
*  The usual "hot model" here..  only with one species for now
*
*/
public class GaussianVLISM extends InterstellarMedium {

	public HelioVector vBulk;
	public double temp, density;
	public static double KB = 1.3807*Math.pow(10,-23);
	public static double MP = 1.672621*Math.pow(10,-27);
	private double Mover2KbT, norm, temp1,temp2,temp3;

	public static double NaN = Double.NaN;

	/**
	* Construct a gaussian VLISM
	*/
	public GaussianVLISM(HelioVector _vBulk, double _density, double _temp) {
		vBulk = _vBulk;
		temp = _temp;
		density = _density;
		//System.out.println("created gaussian vlism with ro: " + density + " temp: " + temp + " bulk: " + vBulk.getR());
		// normalization for maxwellian distribution

		//  trying power of 1/2 to see what the difference is
		norm = density*Math.pow(4.0*MP/2.0/Math.PI/KB/temp,3.0/2.0);
		//norm = density*Math.sqrt(4.0*MP/2.0/Math.PI/KB/temp);
		// argument of exponential requires this factor
		Mover2KbT = 4.0*MP/2.0/KB/temp;
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
		temp1 = v1-vBulk.getX(); temp2 = v2-vBulk.getY(); temp3 = v3-vBulk.getZ();
		temp1 = temp1*temp1+temp2*temp2+temp3*temp3;
		return norm*Math.exp(-Mover2KbT*temp1);
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
	*
	*   redoing testing w/ monte carlo integration - test MultipleIntegration here now 12/2010
	*/
	public static final void main(String[] args) {

		final GaussianVLISM gvli = new GaussianVLISM(new HelioVector(HelioVector.CARTESIAN,
				-27000.0,0.0,0.0), 100.0, 6000.0);
		System.out.println("Created GVLISM" );


/*
		int res = 100;
		file f = new file ("gauss_samp.txt");
		f.initWrite(false);
		for (double v = -40000; v<10000; v+= 50000/res) {
			double dens = gvli.heliumDist(v,0,0);
			f.write(v+"\t"+dens+"\n");
		}
		f.closeWrite();
*/


		//System.out.println("Maxwellian VLISM created, norm = " + norm);
		MultipleIntegration mi = new MultipleIntegration();
		//mi.setEpsilon(0.001);
		FunctionIII dist = new FunctionIII () {
			public double function3(double x, double y, double z) {
				return gvli.heliumDist(x, y, z);
			}
		};

		System.out.println("Test dist: " + dist.function3(20000,0,0));

		//double[] limits = new double[6];
		//limits[0]=mi.MINUS_INFINITY; limits[1]=mi.PLUS_INFINITY;
		//limits[2]=mi.MINUS_INFINITY; limits[3]=mi.PLUS_INFINITY;
		//limits[4]=mi.MINUS_INFINITY; limits[5]=mi.PLUS_INFINITY;
		double[] limits2 = new double[6];
		limits2[0]=-100000.0; limits2[1]=100000.0;
		limits2[2]=-100000.0; limits2[3]=100000.0;
		limits2[4]=-100000.0; limits2[5]=100000.0;
		//double z = mi.integrate(dist,limits2,1000);
		//o(" plain integration: " + z);
		double z2 = mi.mcIntegrate(dist,limits2,100000);
		o(" mc at 100000: " + z2);
		double z3 = mi.mcIntegrate(dist,limits2,1000000);
		o(" mc at 1000000: " + z3);
		//System.out.println("answer:  plain, mc1, mc2 = " + z+ " "+ z2 + "  " + z3);


		FunctionIII distZw = new FunctionIII () {
			public double function3(double r, double p, double t) {
				HelioVector hv = new HelioVector(HelioVector.SPHERICAL,r,p,t);
				return gvli.heliumDist(hv)*r*r*Math.sin(t);
			}
		};
		double[] limitsZw = {0.0,100000.0,0.0,2*Math.PI,0,Math.PI};
		double ans = mi.mcIntegrate(distZw, limitsZw,1000000);

		System.out.println("answer in spherical coords = " + ans );


	}

	public static void o(String s) {
		System.out.println(s);
	}

}
