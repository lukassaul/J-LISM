
import cern.jet.math.Bessel;
import java.lang.Math;
import java.util.Date;
import drasys.or.nonlinear.*;

public class BofTheta {
	public static double Ms = 2 * Math.pow(10,30);
	public static double G = 667 * Math.pow(10,-13);
	public static double AU = 15* Math.pow(10,12);

	private double r, v0;
	private Simpsons s;
	private EquationSolution es;

	public BofTheta(double r, double v0) {
		this.r = r; this.v0 = v0;
		System.out.println("r= " + r);
		s = new Simpsons(); s.setEpsilon(.001);
		es = new EquationSolution(); es.setAccuracy(1); es.setMaxIterations(100);
		System.out.println(es.getAccuracy() + " is accurate?");
		System.out.println("r= " + r);
		System.out.println(es.getMaxIterations()+"");
	}

	class DThetaDr implements FunctionI {
		private double b = 0;
		public DThetaDr(double b) {
			this.b = b;
		}
		public double function(double rV) {
			//System.out.println("Called function- "+b+ " " + r + " " +rV);
			return v0*b/(rV*rV)*Math.sqrt(2*Ms*G/rV - v0*v0*b*b/(rV*rV) + v0*v0);
		}
	}

	class ThetaOfB implements FunctionI {
		public double function(double b) {
			//System.out.println("Calling ThetaOfB");
			DThetaDr dthDr = new DThetaDr(b);
			//System.out.println("trying integral..." + dthDr.b);
			try {
				double q = s.integrate(dthDr, 1000*AU, r);
				System.out.println("Theta= " + q);
				return q;
			}catch (Exception e) {
				e.printStackTrace();
				return 0;
			}
		}
	}


	public double calculate(double theta) {
		ThetaOfB tob = new ThetaOfB();
		try {
			double bLimit = Math.sqrt(2*Ms*G*r/(v0*v0) + r*r);
			return es.solve(tob, -bLimit, bLimit, theta);
		}catch (Exception e) {
			e.printStackTrace();
			return 0;
		}
	}
}

