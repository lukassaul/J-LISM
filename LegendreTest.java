import java.util.Date;
public class LegendreTest {

	public static void main(String[] args) {

		Date d1 = new Date();
		System.out.println(Legendre.compute(3,3));
		Date d2 = new Date();
		System.out.println(d2.getTime() - d1.getTime());

		Legendre l = new Legendre();


		file f = new file("deltadata6.txt");
		f.initWrite(false);
		//double cdels[] = new double[100];
		//Arrays.fill(cdels,0.0);
		for (double j=-1.0; j<1.0; j+=+0.04) {
			double delt=0.0;
			for (int i=0; i<50; i++) {
				delt += (2*i+1)/2*l.p(i,j)*l.p(i,0.0);
			}
			//System.out.println("testing");
			System.out.println(j+"\t"+delt);
			f.write(j+"\t"+delt+"\n");
		}
		f.closeWrite();
		System.exit(0);

	}
}
