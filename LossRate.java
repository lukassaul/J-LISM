import java.io.File;
import java.lang.Math;
import java.util.*;


public class NeutralDistribution {

	public static double Ms = 2 * 10^30;
	public static double G = 6.67 * .00000000001;
	public static double K = 1; // correct this soon!
	public static double Mhe = 1 * 10^-80;

	public double Q, F1, F2, Temp, Vb, N, Vtemp, V0;

	// V0 = Moebius/Rucinski's Vinfinity
	// lowercase parameters are inputs to distribution function

	public NeutralDistribution(double bulkVelocity, double temperature, double density) {
		Temp = temperature;
		Vb = bulkVelocity;
		N = density;
		Vtemp = Math.sqrt(2*K*Temp/Mhe);  //average velocity due to heat
	}

	//NOTE:  for now parameters of distribution must be given in coordinates
	//    relative to inflow direction

	public double distributionFunction(double r, double theta,
									   double vr, double vt, double vPhi) {
		// Damby & Camm 1957
		Q = Math.sqrt(G*Ms/r);
		V0 = Math.sqrt(vr*vr + vt*vt - 2*Q*Q);
		F1 = 2*Vb*vr/(Vtemp*Vtemp) * Math.cos(theta)*(V0 - vr)*(V0 - vr) / (V0*(V0-vr)+Q*Q) -
				((Vb*Vb)+(V0*V0)+2*V0*Vb*Math.cos(theta))/(Vtemp*Vtemp);
		F2 = 2*V0*Vb/(vt*vt) * Math.sin(theta)*(V0-vr)*vt/ (V0*(V0-vr) + (Q*Q));
		return N*Math.pow((Math.sqrt(Math.PI)*Vtemp),-3) * Math.exp(F1 + F2*Math.sin(vPhi));
	}

}