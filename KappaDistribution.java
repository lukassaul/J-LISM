/**
*  A basic power law in ENERGY
*
*/
public class PowerLawVLISM extends InterstellarMedium {

	public HelioVector vBulk;
	public double power, density, lim1, lim2, norm, ee;
	public static double NaN = Double.NaN;
	public static double MP = 1.672621*Math.pow(10,-27);
	public static double EV = 6.24150965*Math.pow(10,18);


	/**
	* Construct a power law VLISM
	*/
	public PowerLawVLISM(HelioVector _vBulk, double _density, double _power, double _ll, double _hh) {
		vBulk = _vBulk;
		power = _power;
		density = _density;
		lim1 = _ll;
		lim2 = _hh;
		// calculate normalization
		norm = density/4/Math.PI/ (Math.pow(lim2,power+3.0)- Math.pow(lim1,power+3.0))*(power+3.0);
		norm=norm;
		System.out.println("norm: " + norm);
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
		// calculate energy
		double vv1 = v1-vBulk.getX();
		double vv2 = v2-vBulk.getY();
		double vv3 = v3-vBulk.getZ();
		double vv = Math.sqrt(vv1*vv1+vv2*vv2+vv3*vv3);
	//	ee = 4*MP/2*(vv1*vv1+vv2*vv2+vv3*vv3)*EV;
	//	System.out.println("ee: " + ee);
		if (vv>lim1 && vv<lim2) return norm*Math.pow(vv,power);
		else return 0.0;
	}

	/**
	* default - pass in 1 energy
	*/
	//public double dist(double energy) {
	//	if (v1==NaN | v2==NaN | v3==NaN) return 0.0;
		// calculate energy
		//double vv1 = v1-vBulk.getX();
		//double vv2 = v2-vBulk.getY();
		//double vv3 = v3-vBulk.getZ();
		//ee = 4*MP/2*(vv1*vv1+vv2*vv2+vv3*vv3)*EV;
	//	System.out.println("ee: " + ee);
	//	ee=energy;
	//	if (ee>lim1 && ee<lim2) return norm*Math.pow(ee,power);
	//	else return 0.0;
	//}

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
	*
	*/
	public static final void main(String[] args) {

		final PowerLawVLISM gvli = new PowerLawVLISM(new HelioVector(HelioVector.CARTESIAN,
				-20000.0,0.0,0.0), 1.0, -4.0, 1000.0, 300000.0);
		System.out.println("Created PLVLISM" );

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
		//double z = mi.mcIntegrate(dist,limits,128  );
		double z2 = mi.mcIntegrate(dist,limits2,1000000);

		System.out.println("Integral should equal density= "  + z2);
	}

}
