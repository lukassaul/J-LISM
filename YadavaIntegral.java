
//import cern.jet.math.Bessel;
import java.lang.Math;
import java.util.Date;

public class YadavaIntegral {

	public YadavaIntegral() {
		MultipleIntegration mi = new MultipleIntegration();

		double[] lims2 = new double[4];
		lims2[0]=0.0;
		lims2[1]=Math.PI;
		lims2[2]=0.0;
		lims2[3]=Math.PI;


		FunctionII yadavaEnergyLoss = new FunctionII () {
			double R = 2.8E-15;
			public double function2(double b, double theta) {
				return  Math.sin(theta) * (1-R*Math.sin(theta)/b) / b / Math.exp(R*Math.sin(theta)/b);
			}
		};

		double test2 = mi.integrate(yadavaEnergyLoss,lims2);

		System.out.println(test2+"");
	}

	public static void main(String[] args) {
		YadavaIntegral yi = new YadavaIntegral();
	}
}

