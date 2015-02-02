
public class TestDistribution {

   public static void main(String[] args) {
	    double AU = 15 * 10^13; // m
	    SimpleNeutralDistribution snd = new SimpleNeutralDistribution(
		   		25000, 0, 0, 10000, 7*Math.pow(10,-3), 0, 6.8*Math.pow(10,-8));

		file f=new file("vrDist.txt");
		f.initWrite(false);
		for (int v0=20000; v0<40000; v0=v0+1000) {
			double q = snd.N(AU, 135*Math.PI/180, v0);
			f.write(v0+" "+q+" " +System.getProperty("line.separator"));
		}
		f.closeWrite();
	}
}
