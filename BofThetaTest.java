import java.util.Date;
public class BofThetaTest {

	public static void main(String[] args) {
		//System.out.println(10^2+" "+10^1+" "+1000+" "+Math.pow(10,2));

		double AU = 1.5* Math.pow(10,11);
		Date d1 = new Date();
		VIonBulk vib = new VIonBulk(AU, 28000,135*Math.PI/180);
		System.out.println(vib.Vr()+" " + vib.Vperp());
		Date d2 = new Date();
		System.out.println("Took: "+(d2.getTime()-d1.getTime())+ " milliseconds");
	}
}
