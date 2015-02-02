
import java.util.StringTokenizer;

public class IbexFitter {

	public CurveFitter cf;
	public IBEXLO_09fit ifit;

	public IbexFitter() {

		file f = new file("first_pass.txt");
		int numLines = 926;
		numLines=numLines/3;
		f.initRead();
		double[] x = new double[numLines];
		double[] y = new double[numLines];
		String line = "";
		for (int i=0; i<numLines; i++) {
			line = f.readLine();
			StringTokenizer st = new StringTokenizer(line);
			x[i]=Double.parseDouble(st.nextToken());
			//x[i]-=1.0; // IMPORTANT _ SHIFT BY 1.0 TO MOVE TO HELIOCENTRIC FRAME!!
			y[i]=Double.parseDouble(st.nextToken());
			//System.out.println("row: " + i + " x: " + x[i] + " y: " + y[i]);
			line = f.readLine();
			line = f.readLine();
		}
		f.closeRead();


		// loop through parameters here
		// standard params:
		double lamda = 74.7;
		double v = 26000.0;
		int res = 2;
		double lamdaWidth = 15;
		double vWidth = 10000;
		double lamdaMin = 65.0;
		double vMin=21000.0;

		double lamdaDelta = lamdaWidth/res;
		double vDelta = vWidth/res;

		lamda = lamdaMin;
		v = vMin;
		file outF = new file("testOut.txt");
		outF.initWrite(false);
		for (int i=0; i<res; i++) {
			for (int j=0; j<res; j++) {


				System.out.println("starting doFit routine with i=" + i + " j= " + j);

				cf = new CurveFitter(x,y,lamda,v);
				cf.doFit(CurveFitter.IBEX1);

				outF.write(lamda+"\t"+v+"\t");
				outF.write(cf.minimum_pub+"\t");
				outF.write(cf.bestParams[0]+"\n");

				// get values of y and z at minimum
				//double[] param = min.getParamValues();
				//bestParams = param;
				v+=vDelta;

			}
			lamda+=lamdaDelta;

		}

		outF.closeWrite();



	}


	public static final void main(String[] args) {
		IbexFitter theMain = new IbexFitter();
	}
}