/**
*
* Use this to calculate the v dist. of int. helium using
*   the assumption of isotropic adiabatic cooling
*
*/
public class AdiabaticModel  {

	public double AU  = 1.49598* Math.pow(10,11); //meters
	public double NaN = Double.NaN;
	public double MAX = Double.MAX_VALUE;
	public double N0  = Math.pow(10,-6);
	public double GAMMA = 3.0/2.0;  // energy dependence in eflux ?? - inverse gamma in derivation..
	//public double beta = 1.0*Math.pow(10,-7); // that is for He
	public double beta = 7.4*Math.pow(10,-7);
	static double PI = Math.PI;

	public double VSW = 350000.0;
	public static double v_infinity = 27500.0;  //m/s
	public double G = 6.673 * Math.pow(10,-13);  // m^3/s^2/kg
	public static double Ms = 1.98892 * Math.pow(10,30);  // kg
	public static double VINF = 26000.0;

    public AdiabaticModel() {

	}

	public double f(double norm, double gam, double w) {
		GAMMA = gam;
		return f(w)*norm;
	}

	public double f(double norm, double gam, double w, double _beta) {
		beta = _beta;
		return f(norm,gam,w);
	}

	public double f(double w) {
		//if (w<=1-Math.sqrt(VINF*VINF+2*Ms*G/AU)/VSW)
			return Math.pow(w,1.0-1.0/GAMMA)*N(AU*Math.pow(w,1.0/GAMMA))/GAMMA;
		//return 0.0;
	}

	public double f_old(double w) {
		if (w<=1) return Math.pow(w,1.0-1.0/GAMMA)*N(AU*Math.pow(w,1.0/GAMMA));
		return 0.0;
	}

	public double f_test(double w) {
		return Math.pow(w,1.0-1.0/GAMMA) * N(AU*Math.pow(w,1.0/GAMMA))
				/ square(1.0 + 1/VSW*Math.sqrt(VINF*VINF+2*Ms*G/AU/Math.pow(w,1.0/GAMMA)));
	}


	public double f(double gam, double w) {
		GAMMA = gam;
		return f(w);
	}

	public static double square(double s) {return s*s;}

	/**
	*
	*Vasyliunas & Siscoe for comparison
	*
	*/
	public double f_vs1(double w) {
		return Math.exp(-Math.pow(w,-2.0));
	}

	public double f_vs2(double w) {
		double lam=15.2;
		return Math.pow(w,-0.75)*Math.exp(-lam*Math.pow(w,-1.5));
	}



	/**
	* The model of interstellar neutral density..
	*
	*   based on cold model, due upwind..
	*   see Moebius SW&CR homework set #4!
	*/
	private double N(double r) {
		return N0*Math.exp(-beta*AU*AU*Math.sqrt(v_infinity*v_infinity+/*2.0**/G*Ms/2.0/r)/G/Ms);
		//return N0*Math.exp(-2*beta*AU*AU/r/v_infinity);
	}

	public static double e(double a, double b) {
		return a*Math.pow(10,b);
	}

	public static final void main(String[] args) {
		AdiabaticModel sm = new AdiabaticModel();
		sm.G=Math.pow(10.0,-11.0);

		file f = new file("adiabatic_out_close2.dat");
		f.initWrite(false);
		for (double w=0.0; w<=3.0; w+=0.02) {
			f.write(w+"\t");
			f.write(sm.f(w)+"\t"+sm.f_vs1(w)+"\t"+sm.f_vs2(w)+"\n");

			//f.write(""+sm.f(1.0,1.5,w));
			//f.write(sm.f(1.0,1.5,w,e(5.0,-9.0))+"\t"+sm.f(1.0,1.5,w,e(1.0,-8.0))+"\t"+sm.f(1.0,1.25,w,e(5.0,-8.0))+"\t"+
			//		sm.f(1.0,1.5,w,e(1.0,-7.0))+"\t"+sm.f(1.0,1.5,w,e(5.0,-7.0))+"\t"+sm.f(1.0,1.25,w,e(1.0,-6.0))
			//				);
			//f.write("\n");
		}
		f.closeWrite();


		/*double[] test = new double[20];
		double[] test1 = new double[20];
		double[] test2 = new double[20];
		double[] test3 = new double[20];
		double[] test4 = new double[20];
		double[] test5 = new double[20];
		double[] test6 = new double[20];
		double[] gamma = {0.3d, 0.8d, 1.3d, 1.8d};
		file f = new file("adiabatic_out_test2.dat");
		f.initWrite(false);
		int index = 0;
		for (double w=0.0; w<=1.0; w+=0.02) {
			index=0;
			f.write(w+"\t");
			while (index<4) {
				f.write(sm.f(gamma[index],w)+"\t");
				index++;
			}
			f.write("\n");
		}
		f.closeWrite();
		*/
	}
}
