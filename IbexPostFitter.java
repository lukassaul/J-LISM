
import java.util.StringTokenizer;

public class IbexPostFitter {

	public CurveFitter cf;

	public IbexPostFitter() {

		file asciiOutFile = new file("second_pass_ascii.txt");

		file f = new file("dirD.txt");
		int numLines = 900;

		// load the real data
		file df = new file("second_pass.txt");
		int numLines2 = 918;
		df.initRead();
		double[] xD = new double[numLines2];
		double[] yD = new double[numLines2];
		String line = "";
		for (int i=0; i<numLines2; i++) {
			line = df.readLine();
			StringTokenizer st = new StringTokenizer(line);
			xD[i]=Double.parseDouble(st.nextToken());
			//x[i]-=1.0; // IMPORTANT _ SHIFT BY 1.0 TO MOVE TO HELIOCENTRIC FRAME!!
			yD[i]=Double.parseDouble(st.nextToken());
			System.out.println("row: " + i + " x: " + xD[i] + " y: " + yD[i]);
			//line = f.readLine();
			//line = f.readLine();
		}
		df.closeRead();

		// then first lets get file list of model data
		String fileList[] = new String[900];
		f.initRead();
		line = "";
		String garbage = "";
		int ii=0;
		while ((line=f.readLine())!=null) {
			StringTokenizer st1 = new StringTokenizer(line);
			// throw out 4 and take the 5th
			garbage = st1.nextToken();
			garbage = st1.nextToken();
			garbage = st1.nextToken();
			garbage = st1.nextToken();
			fileList[ii]=st1.nextToken();
			//System.out.println(i+" " + fileList[i]);
			ii++;
		}

		// now get the params from the model filename
		double[] x = new double[30];
		double[] y = new double[30];
		double[] z = new double[900];
		int xNum = 0;
		int yNum = 0;
		for (int j=0; j<fileList.length; j++) {
			StringTokenizer st2 = new StringTokenizer(fileList[j],"_vl");
			garbage = st2.nextToken(); o(garbage);
			String sss = st2.nextToken(); o("good?: " + sss);
			double test1 = Double.parseDouble(sss);
			//System.out.println(test1);
			if (!contains(x,test1)) {
				x[xNum]=test1;
				xNum++;
			}


			String aNumber = st2.nextToken();
			o("good2?: " + aNumber);
			int len = aNumber.length();
			aNumber = aNumber.substring(0,len-4);
			double test2 = Double.parseDouble(aNumber);
			if (!contains(y,test2)) {
				y[yNum]=test2;
				yNum++;
			}
			//System.out.println(xNum+" "+yNum);

			// ok we have params now lets load the file and get the fit factor
			cf = new CurveFitter(xD,yD,fileList[j]);
			cf.doFit(CurveFitter.IBEXFIT);
			z[j] = (cf.minimum_pub);
			System.out.println("did fit: "+ j + " : " + z[j]);


		}

		// write output to ascii file for later graphing
		asciiOutFile.initWrite(false);
		for (int i=0; i<x.length; i++) 	asciiOutFile.write(x[i]+"\n");
		for (int i=0; i<y.length; i++) 	asciiOutFile.write(y[i]+"\n");
		for (int i=0; i<z.length; i++) 	asciiOutFile.write(z[i]+"\n");
		asciiOutFile.closeWrite();

		JColorGraph jcg = new JColorGraph(x,y,z);

		String unitString = "log (sum square model error)";
		jcg.setLabels("IBEX-LO","2010",unitString);

		jcg.run();
		jcg.showIt();





	}

	public static boolean contains(double[] set, double test) {
		for (int i=0; i<set.length; i++) {
			if (set[i]==test) return true;
		}
		return false;
	}

	public static void o(String s) {
		System.out.println(s);
	}

	public static final void main(String[] args) {
		IbexPostFitter theMain = new IbexPostFitter();
	}
}