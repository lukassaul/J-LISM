
import cern.jet.math.Bessel;
import java.lang.Math;
import java.util.Date;
import drasys.or.nonlinear.*;


/**  This class computes the pickup distribution as a function of
*    vr, vt, ( and r, theta)
*/

public class 2DVelocityDistribution {
	// constants not correct yet
	public static double Ms = 2 * Math.pow(10,30);
	public static double G = 6.67 * Math.pow(10,-11);
	public static double AU = 1.5* Math.pow(10,11);
	public static double K = 1.380622 * Math.pow(10,-23);    //1.380622E-23 kg/m/s/s/K
	public static double Mhe = 4*1.6726231 * Math.pow(10,-27); // 4 * 1.6726231 x 10-24 gm

	public double Q, F1, F2, Temp,  Vb, N, Vtemp, V0, mu, r, theta, betaPlus, betaMinus;
	public double phi0, theta0;  // direction of bulk flow - used when converting to GSE

	//public Trapezoidal s;
	public Integration s;
	public NeutralDistribution nd; // this is our neutral helium distribution
	public SurvivalProbability sp; // % of neutral ions still left at r,v
	public 2DVelocityDistribution(double bulkVelocity, double phi0, double theta0,
						double temperature, double density,
						double radiationToGravityRatio, double lossRate1AU) {
		// the constructor just sets up the parameters
		Temp = temperature;
		Vb = bulkVelocity;
		this.phi0 = phi0;
		this.theta0 = theta0;
		N = density;
		mu = radiationToGravityRatio;
		betaMinus = lossRate1AU;
		Vtemp = Math.sqrt(2*K*Temp/Mhe);  //average velocity due to heat

		sp = new SurvivalProbability(betaMinus);
		nd = new NeutralDistribution(Vb, phi0, theta0,
						Temp, N, mu);
	}

	public double N(double r, double theta, double vr, double vt) {

		V0 = Math.sqrt(vr*vr + vt*vt - 2*Q*Q);
		F1 = 2*Vb*vr/(Vtemp*Vtemp) *
		           Math.cos(theta)*(V0 - vr)*(V0 - vr) / (V0*(V0-vr)+Q*Q) -
			           ((Vb*Vb)+(V0*V0)+2*V0*Vb*Math.cos(theta))/(Vtemp*Vtemp);
		F2 = 2*V0*Vb/(Vtemp*Vtemp) * Math.sin(theta)*(V0-vr)*vt/ (V0*(V0-vr) + (Q*Q));

		return sp.compute(r,theta,vt,vr) * N*Math.pow((Math.sqrt(Math.PI)*Vtemp),-3) *
				Math.exp(F1)*Bessel.i0(F2);

	}

	// this function returns N(r,theta, vt)
	public double N(double r, double theta, double vt) {
		this.r=r; this.theta=theta;

		Q = Math.sqrt((1-mu)*G*Ms/r);
		Nvr nvr = new Nvr(vt);
		//s = new Trapezoidal(); s.setMaxIterations(15); s.setEpsilon(1000000);
		s=new Integration();
		s.setFunction(nvr);
		s.setMaxError(.01);

		try {
			return s.integrate( -100000, 100000);
		}catch(Exception e) {
			e.printStackTrace();
			return 0;
		}
	}


	class Nvr implements FunctionI, Function {
		// this is the function we integrate over vr
		public double vt;
		public Nvr(double vt_) {
			System.out.println("nvr constructor");
			vt = vt_;
	//		System.out.println("set vt");
		}
		public double function(double vr) {
	//		System.out.println("Trying nvr.function vr= " + vr);
			V0 = Math.sqrt(vr*vr + vt*vt - 2*Q*Q);
			F1 = 2*Vb*vr/(Vtemp*Vtemp) * Math.cos(theta)*(V0 - vr)*(V0 - vr) / (V0*(V0-vr)+Q*Q) -
								((Vb*Vb)+(V0*V0)+2*V0*Vb*Math.cos(theta))/(Vtemp*Vtemp);
			F2 = 2*V0*Vb/(Vtemp*Vtemp) * Math.sin(theta)*(V0-vr)*vt/ (V0*(V0-vr) + (Q*Q));

			return sp.compute(r,theta,vt,vr)*N*Math.pow((Math.sqrt(Math.PI)*Vtemp),-3) *
					Math.exp(F1)*Bessel.i0(F2);

			// here we integrate over L(r,v)*N(r,v) with L = survival probability.
			// this could be all put into this routine but I leave it out for other
			// factors to be added in later
		}
	}
}