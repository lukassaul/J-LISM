import java.util.Date;
import java.lang.Math;
public class NTest {

	public static void main(String[] args) {
		Date d1 = new Date();
		//NeutralDistribution nd = new NeutralDistribution( 12, 123, 142);

		double AU = 15 * 10^13; // m
		SimpleNeutralDistribution snd = new SimpleNeutralDistribution(
				25000, 0, 0, 10000, 7*Math.pow(10,-6), 0, 6.8*Math.pow(10,-8));



			double q = snd.N(AU, 135*Math.PI/180, 25000);

			System.out.println(q+"");



		/*double[] array;
		array = new double[100000];
		for (int i = 1; i<100000; i++) {
			array[i] = i;
		//	double test = nd.distributionFunction(123123, .21, 100+i, 421, .13);
		}*/
		Date d2 = new Date();
		System.out.println(d2.getTime() - d1.getTime()+"");
	}
}


