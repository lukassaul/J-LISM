import java.util.Date;
public class ThetaOfBTest {

	public static void main(String[] args) {
		//System.out.println(10^2+" "+10^1+" "+1000+" "+Math.pow(10,2));

		double AU = 1.5* Math.pow(10,11);
		Date d1 = new Date();
		BofTheta vib = new BofTheta(AU, 28000);
		System.out.println(vib.testThetaOfB(AU));
		System.out.println(vib.testThetaOfB(3*AU));
		System.out.println(vib.testThetaOfB(1.8*AU));
		System.out.println(vib.testThetaOfB(2*AU));
		Date d2 = new Date();
		System.out.println("Took: "+(d2.getTime()-d1.getTime())+ " milliseconds");
	}
}
