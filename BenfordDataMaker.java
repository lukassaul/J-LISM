import java.util.Random;
public class BenfordDataMaker {
	public file f;
	public Random r;
	public BenfordDataMaker() {
		r = new Random();
    	f = new file("testdatarandom.txt");
    	f.initWrite(false);
    	for (int i =0; i<10000; i++) {
			f.write(r.nextDouble()+"\t");
			f.write(r.nextGaussian()+"\t");
			f.write(nextLogNormal()+"\n");
		}
		f.closeWrite();
	}

	public double nextLogNormal() {
		double q = r.nextGaussian();
		return Math.exp(q);
	}

	public static final void main(String[] args) {
		BenfordDataMaker bdm = new BenfordDataMaker();
	}
}