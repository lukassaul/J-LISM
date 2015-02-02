
/**
* Use a NeutralDistribution to compute a NeutralDensity in the heliosphere
*
*  this is done by integrating over all velocity space for a desired grid
*
*  x, y, built around origin from + to - 'scale'.
*/
public class NeutralDensity {
	public int mc_num = 80000;
	private static int gridsize = 40; // gridsize^2 spatial calcualtion points
	public static double AU = 1.49598* Math.pow(10,11); //meters
	private static double scale = 100.0; // edge size in AU
	public String filename = "kappa_2_20000_40.dat";

	//private static double vMax = 10000000.0 // 10 k km/s


	// x and y are now the same but we hang onto them for clarity
	private double[] x = new double[gridsize];
	private double[] y = new double[gridsize];
	private double[] z = new double[gridsize*gridsize];

	private NeutralDistribution nd;
	private MultipleIntegration mi;


	/**
	*  Default, add lism params by hand here
	*/
	public NeutralDensity () {

/*		this(new NeutralDistribution(
			new KappaVLISM(
				new HelioVector(HelioVector.CARTESIAN,26000.0,0.1,0.1),
				100.0,8000.0,1.6),
			0.0, 0.0));//7*Math.pow(10.0,-8.0)));
*/		
		this(new NeutralDistribution(
			new GaussianVLISM(
				new HelioVector(HelioVector.CARTESIAN,26000.0,0.1,0.1),
				0.015*100.0*100.0*100.0,12500),
			0.0, 0.0));
			
//		
//	public GaussianVLISM gv1;
//bulkHE = new HelioVector(HelioVector.SPHERICAL, 26000.0, (74.68)*Math.PI/180.0, 95.5*Math.PI/180.0);
		//HelioVector bulkO1 = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
//		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,0.015*100.0*100.0*100.0,6000.0);

	}


	/**
	*  Build this object from a distribution.
	*  0th order moment creation here
	*/
	public NeutralDensity(NeutralDistribution nd_) {
		nd = nd_;
		mi = new MultipleIntegration();
		//mi.setEpsilon(0.0001);
		nd.debug=false;

		//calculateGrid();
	}

	private void calculateGrid() {

		double delta = scale*AU/(double)gridsize;
		System.out.println("delta: " + delta);
		System.out.println("gridsize " + gridsize);
		int index = 0;

		for (int i=1; i<=gridsize; i++) {

			x[index]=scale*AU*0.05+delta*i;
			//x[index]=0.0-(scale*AU/2)+delta*i;
			y[index]=0.0-(scale*AU/2)+delta*i;
			System.out.println("x["+index+"] : " + x[index]);
			index++;
		}

		index = 0;
		for (int i=0; i<x.length; i++) {
			for (int j=0; j<y.length; j++) {
				HelioVector pos = new HelioVector(HelioVector.CARTESIAN, x[i], y[j], 0.0);
				z[index]=vint(pos);
				index++;
				System.out.println("calculated z[" + (index-1) + "] : " + z[index-1]);
			}
		}

		file f = new file(filename);
		f.initWrite(false);
		for (int i=0; i<x.length; i++) {
			f.write(x[i]+"\n");
		}
		f.write("\n");
		for (int i=0; i<y.length; i++) {
			f.write(y[i]+"\n");
		}
		f.write("\n");
		for (int i=0; i<z.length; i++) {
			f.write(z[i]+"\n");
		}
		f.closeWrite();


		//JColorGraph jcg = new JColorGraph(x,y,z);
		//jcg.run();
	}


	/**
	* Do the integration!  Pass in a the position as a heliovector
	*/
	private double vint(HelioVector hp) {
		final HelioVector hpf = hp;
		FunctionIII dist = new FunctionIII () {
			// r,p,t here actually a velocity passed in...
			public double function3(double r, double p, double t) {
				return r*r*Math.sin(t)*nd.dist(hpf,new HelioVector(HelioVector.SPHERICAL, r, p, t));
			}
		};
		double[] limits = new double[6];
		limits[0]=0.0; limits[1]=100000;
		limits[2]=0.001; limits[3]=2*Math.PI;
		limits[4]=0.001; limits[5]=Math.PI;
		//limits[4]=Math.sqrt(2*nd.gmodMs/hp.getR()); limits[5]=300000;
		return mi.mcIntegrateSS(dist,limits, mc_num);
	}

	/**
	* Do the integration!  Pass in a the position as a heliovector and a v_r
	*/
	private double vrint(HelioVector hp, double v__r) {
		final HelioVector hpf = hp;
		final double v_r = v__r;
		FunctionII dist = new FunctionII () {
			// r,p,t here actually a velocity passed in...
			public double function3(double p, double t) {
				return v_r*v_r*Math.sin(t)*nd.dist(hpf,new HelioVector(HelioVector.SPHERICAL, v_r, p, t));
			}
		};
		double[] limits = new double[6];
		//limits[0]=0.0; limits[1]=100000;
		limits[0]=0.001; limits[1]=2*Math.PI;
		limits[2]=0.001; limits[3]=Math.PI;
		//limits[4]=Math.sqrt(2*nd.gmodMs/hp.getR()); limits[5]=300000;
		return mi.mcIntegrate(dist, limits, mc_num);
	}

	/**
	* Use for testing first, then for starting routine to generate grid
	*/
	public static void main(String[] args) {
		// first we do some tests		

		// first test - density at 1000AU,1000AU,1000AU

		/*HelioVector tv = new HelioVector(HelioVector.CARTESIAN, 100*AU, 100*AU, 100*AU);
		NeutralDistribution ndist = new NeutralDistribution(
						new GaussianVLISM(new HelioVector(HelioVector.CARTESIAN,50.0,0.0,0.0),
	        					  50000,8000.0)
	        		,0.0,0.0);

		NeutralDensity nd = new NeutralDensity(ndist);

	   	 System.ouct.println("test: " + nd.vint(tv));
		System.out.println("nd.counter: " + ndist.counter);  */

		//HelioVector bulk = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		//HelioVector bulk = new HelioVector(HelioVector.CARTESIAN, 23500.0, 0.0, 0.0);
		//NeutralDistribution ndH = new NeutralDistribution(new GaussianVLISM(bulk,1.0,6500.0),  0.0, 1.0*Math.pow(10,-7));
		//NeutralDistribution ndH = new NeutralDistribution(new KappaVLISM(bulk,1.0,6500.0,10.0),  0.0, 1.0*Math.pow(10,-7));

		//ndH.debug=false;
		//NeutralDensity nd = new NeutralDensity();

		/*NeutralDensity nd = new NeutralDensity(new NeutralDistribution(
			new KappaVLISM(
				new HelioVector(HelioVector.CARTESIAN,26000.0,0.1,0.1),
				0.015*100.0*100.0*100.0,6500,2.0),
			0.0, 0.0));
		nd.filename="kappa_1_6_x5.dat";
		nd.calculateGrid();*/

		NeutralDensity nd = new NeutralDensity(new NeutralDistribution(
			new KappaVLISM(
				new HelioVector(HelioVector.CARTESIAN,26000.0,0.1,0.1),
				0.015*100.0*100.0*100.0,6500,2.1),
			0.0, 0.0));
		nd.filename = "kappa_2_1_x5.dat";		
		//nd.calculateGrid();
		
		double vv=0.0;
		HelioVector hp = new HelioVector(HelioVector.CARTESIAN, -100*AU, 0.01, 0.01);
		file f = new file("veloc_out.dat");
		f.initWrite(false);
 		for (int i=0; i<50; i++) {
			f.write(vv+"\t");
			f.write(nd.vrint(hp,vv)+"\n");
			vv+=2000.0;
		}
		f.closeWrite();

	}
}
