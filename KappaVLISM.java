
import flanagan.analysis.Stat;

/**
*  A basic power law in ENERGY
*
*/
public class KappaVLISM extends InterstellarMedium {

	public HelioVector vBulk;
	public double power, density, kappa, norm, ee, temp, omega;
	public static double NaN = Double.NaN;
	public static double MP = 1.672621*Math.pow(10,-27);
	public static double EV = 6.24150965*Math.pow(10,18);
	public static double KB = 1.3807*Math.pow(10,-23);
	public static double AU = 1.49598* Math.pow(10,11); //meters


	/**
	* Construct a power law VLISM
	*/
	public KappaVLISM(HelioVector _vBulk, double _density, double _temp, double _kappa) {
		vBulk = _vBulk;
		kappa = _kappa;
		density = _density;
		temp = _temp;

		// calculate normalization
		omega = Math.sqrt(2*temp*KB*(kappa-1.5)/kappa/4/MP);
		System.out.println("kappa " + kappa + " omega: " + omega);

		norm = density*gamma(kappa+1.0)/omega/omega/omega/Math.pow(kappa*Math.PI,1.5)/gamma(kappa-0.5);
		System.out.println("norm: " + norm);
	}

	public double gamma(double q) {
		return Stat.gammaFunction(q);
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
		//if (vv>lim1 && vv<lim2) return norm*Math.pow(vv,power);
		//else return 0.0;
		return norm*Math.pow(1+vv*vv/kappa/omega/omega,-kappa-1.0);
	}

	public double dist(HelioVector hv) {
		return dist(hv.getX(), hv.getY(), hv.getZ());
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

		//	public KappaVLISM(HelioVector _vBulk, double _density, double _temp, double _kappa) {


			// ********   Test the kappa function itself *******
			final KappaVLISM gvli = new KappaVLISM(new HelioVector(HelioVector.CARTESIAN,
					-23500.0,0.0,0.0), 1.0, 6500.0, 1.51);
			System.out.println("Created KVLISM" );

			final NeutralDistribution ndk = new NeutralDistribution(gvli, 0.0, 1.0*Math.pow(10,-7));
			ndk.debug=false;


			/** test the distribution at inf and compare to GAUSS  **/
			file f = new file("dist_02au_kap1_5.txt");
			f.initWrite(false);

			final GaussianVLISM ggvli = new GaussianVLISM(new HelioVector(HelioVector.CARTESIAN,
					-23500.0,0.0,0.0), 1.0, 6500.0);


			final NeutralDistribution ndg = new NeutralDistribution(ggvli, 0.0, 1.0*Math.pow(10,-7));
			ndg.debug=false;

			// position
			final HelioVector hpf = new HelioVector(HelioVector.CARTESIAN, 0.2*AU,0.01*AU,0);

			MultipleIntegration mi = new MultipleIntegration();

			double max = 200000.0;
			double logmax = Math.log(max);
			double logmin = Math.log(75000);
			double logstep = (logmax-logmin)/100.0;
			for (double lv=logmin; lv<logmax; lv+= logstep) {
				FunctionIII dist = new FunctionIII () {
					// r,p,t here actually a velocity passed in...
					public double function3(double r, double p, double t) {
						return r*r*Math.sin(t)*ndk.dist(hpf, new HelioVector(HelioVector.SPHERICAL, r, p, t));
					}
				};
				FunctionIII dist2 = new FunctionIII () {
					// r,p,t here actually a velocity passed in...
					public double function3(double r, double p, double t) {
						return r*r*Math.sin(t)*ndg.dist(hpf, new HelioVector(HelioVector.SPHERICAL, r, p, t));
					}
				};
				double[] limits = new double[6];
				limits[0]=Math.exp(lv); limits[1]=Math.exp(lv+logstep);
				limits[2]=0.001; limits[3]=2*Math.PI;
				limits[4]=0.001; limits[5]=Math.PI;
				//limits[4]=Math.sqrt(2*nd.gmodMs/hp.getR()); limits[5]=300000;
				f.write(Math.exp(lv)+"\t"+mi.mcIntegrateSS(dist,limits, 1000000)+"\t"+mi.mcIntegrateSS(dist2,limits, 500000)+"\n");
				System.out.println("did lv: "+ lv);

			}
			f.closeWrite();



		//  ****  Test the Integrals ****
		//System.out.println("Maxwellian VLISM created, norm = " + norm);
		//MultipleIntegration mi = new MultipleIntegration();
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
