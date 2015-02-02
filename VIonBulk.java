
import cern.jet.math.Bessel;
import java.lang.Math;
import java.util.Date;
import drasys.or.nonlinear.*;

public class VIonBulk{
	public static double Ms = 2 * Math.pow(10,30);
	public static double G = 6.67 * Math.pow(10,-11);
	BofTheta boftheta;
	double b,r,v0;

	public VIonBulk(double r, double v0, double theta) {
		this.r = r; this.v0 = v0;
		boftheta = new BofTheta(r, v0);
		b = boftheta.calculate(theta);
	}

	public double Vr() {
		return -Math.sqrt(2*Ms*G/r - v0*v0*b*b/r/r + v0*v0);
	}

	public double Vperp() {
		return v0*b/r;
	}
}




class BofTheta {
	public static double Ms = 2 * Math.pow(10,30);  //kg
	public static double G = 6.67 * Math.pow(10,-11);  // m^3 / (kg s^2)
	public static double AU = 1.5* Math.pow(10,11); // m
	private double r, v0;
	private Trapezoidal s;
	private EquationSolution es;

	public BofTheta(double r, double v0) {
		this.r = r; this.v0 = v0;
		System.out.println("r= " + r);
		s = new Trapezoidal(); s.setEpsilon(.1);
		es = new EquationSolution(); es.setAccuracy(100); es.setMaxIterations(100);
	}

	class DThetaDr implements FunctionI {
		private double b = 0;
		public DThetaDr(double b) {
			this.b = b;
		}
		public double function(double r_) {
			return v0*b/(r_*r_) / Math.sqrt(2*Ms*G/r_ - v0*v0*b*b/(r_*r_) + v0*v0);
		}
	}

	class ThetaOfB implements FunctionI {
		public double function(double b) {
			System.out.println("ThetaOfB called: b="+b);
			DThetaDr dthDr = new DThetaDr(b);
			try {
				double test= s.integrate(dthDr, 1000*AU, r);
				System.out.println("integral returns: " + test);
				return test;
			}catch (Exception e) {
				e.printStackTrace();
				return 0;
			}
		}
	}

	public double testThetaOfB(double b_) {
		ThetaOfB tobTest = new ThetaOfB();
		return tobTest.function(b_);
	}

	public double calculate(double theta) {
		ThetaOfB tob = new ThetaOfB();
		try {
			double bLimit = Math.sqrt(2*Ms*G*r/(v0*v0) + r*r);
			System.out.println("blimit = "+bLimit);
			double q= es.solve(tob, -bLimit+1000, bLimit-1000, theta);
			System.out.println("b (AU)= " + q/AU);
			return q;
		}catch (Exception e) {
			e.printStackTrace();
			return 0;
		}
	}
}

