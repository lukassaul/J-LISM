/**
* A "hot" neutral helium model of the LISM
*   (gaussian)
*/
public class SimpleHeLISM extends InterstellarMedium {

	public static double K = 1.380622 * Math.pow(10,-23);    //1.380622E-23 kg/m/s/s/K
	public static double Mhe = 4*1.6726231 * Math.pow(10,-27); // 4 * 1.6726231 x 10-24 gm

	private double temp, vb, n, vtemp, d;
	private double phi0, theta0;  // direction of bulk flow - used when converting to GSE
	private double norm;
	private double factor;
	/**
	* Default parameters for default constructor
	*
	*/
	public SimpleHeLISM() {
		this(28000.0, 0.0, 0.0, 2000.0, 1.5E-7);
	}

	/**
	* Construct an interstellar neutral helium gaussian distribution
	*parameters:
	*    bulk velocity of VLISM
	*    phi0, phi angle of bulk inflow
	*    theta0, theta angle of bulk inflow
	*    temperature, T(VLISM)
	*	 density (VLISM)
	*    radiationToGravityRation (mu)
	*    ionization at 1AU (beta0)
	*/
	public SimpleHeLISM (double bulkVelocity, double phi0, double theta0,
						double temperature, double density)  {
		vb = bulkVelocity;
		this.phi0 = phi0;
		this.theta0 = theta0;
		temp = temperature;
		d = density;
		factor = -Mhe/2/K/temp;

		//average velocity due to heat
		vtemp = Math.sqrt(-1/factor);

		// don't want to compute this every time...
		norm = d*Math.pow(vtemp*Math.sqrt(Math.PI),-3);
		/*
		// DEBUG HERE
		file f = new file("test_dist.txt");
		f.initWrite(false);
		double minRange = -100000;
		double maxRange = 100000;
		int steps = 100;
		double delta = (maxRange-minRange)/steps;
		for (int i=0; i<steps; i++) {
			f.write((minRange+i*delta)+"\t"+heliumDist((minRange+i*delta),0,0)+"\t");
			f.write(heliumDist(0,(minRange+i*delta),0)+"\t");
			f.write(heliumDist(0,0,(minRange+i*delta))+"\n");
		}
		f.closeWrite();
		*/
	}

	/**
	* Returns the phase space density of neutral helium at this velocity
	*/
	public double heliumDist(double vx, double vy, double vz) {
		double r = Math.sqrt((vx-vb)*(vx-vb) + vy*vy + vz*vz);
		return norm*Math.exp(factor*r*r);
	}


	/**
	* For testing
	*/
	public static final void main(String[] args) {

		SimpleHeLISM shl = new SimpleHeLISM(27000,0,0,8000,1);
	}

}


/*

		double temp = 2000;
		double density = 1;
		double vbulk = 27000;
		double norm = Math.sqrt((4*MP/2/PI/KB/temp)*(4*MP/2/PI/KB/temp)*(4*MP/2/PI/KB/temp));
		return norm*Math.exp(-4*MP*((v1-vbulk)*(v1-vbulk)+v2*v2+v3*v3)/2/KB/temp);

		*/




