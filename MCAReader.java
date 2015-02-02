
import java.io.*;
/**
* Read and parse a single "*.MCA" file
*  These MCA files are outputs of our positional MCP detector
*  sometimes known as "quanta" or "carlson" detector
*
*/
public class 3DFileReader {

	public DataInputStream dis;


	int maxbin = 128;
	int betamax = 3;
	int alphamax = 109;
	int alphabegin = -7;
    float pictratio = alphamax/betamax;
	float[][] image = new float[maxbin][maxbin];


	public 3DFileReader (String t) {
		String testFile = t;
		try {
			dis = new DataInputStream(new FileInputStream(testFile));
			short alpha = dis.readShort();
			float beta = dis.readFloat();
			float theta = dis.readFloat();
			float bm = dis.readFloat();
			float energy = dis.readFloat();
			float acctime = dis.readFloat();
			float totalCounts = dis.readFloat();

			// end reading header
			for (int i=0; i<maxbin; i++) {
				for (int j=0; j<maxbin; j++) {
					image[i][j]=dis.readFloat();
				}
			}

			System.out.println("read file: " + testFile);
			System.out.println("alpha: " + alpha);
			System.out.println("beta: " + beta);
			System.out.println("theta: " + theta);
			System.out.println("bm: " + bm);
			System.out.println("acctime: " + acctime);
			System.out.println("total counts: " + totalCounts);
			System.out.println("0,0 : " + image[0][0]);

		}
		catch (Exception e) {
			e.printStackTrace();
			System.out.println("Couldn't find file: " + testFile);
		}

	}

	public static final void main(String[] args) {
		if (args.length<1) {3DFileReader mr = new 3DFileReader("1209_1.mca");}
		else {3DFileReader mr = new 3DFileReader(args[0]);}
	}

}