
 public class IbexParamPlotter {
 	public CurveFitter cf;
 	public double[] x,y;  // the data we are fitting to
 	public IbexParamPlotter() {

 	}


 	/**
    *  Use this to fit a curve or
    * for testing...
    *
    *  Reads a file for data and fits to a curve at command line
    *
    */
    public static final void main(String[] args) {

		int numLines = 566;
		int linesToSkip = 0;



		file f = new file("spring_data_filtered.txt");
		f.initRead();
		double[] x = new double[numLines];
		double[] y = new double[numLines];
		String line = "";

		for (int i=0; i<linesToSkip; i++) {
			line = f.readLine();
		}

		for (int i=0; i<numLines; i++) {
			line = f.readLine();
			StringTokenizer st = new StringTokenizer(line);
			x[i]=Double.parseDouble(st.nextToken());
			//x[i]-=1.0; // IMPORTANT _ SHIFT BY 1.0 TO MOVE TO HELIOCENTRIC FRAME!!
			y[i]=Double.parseDouble(st.nextToken());
			//System.out.println("row: " + i + " x: " + x[i] + " y: " + y[i]);
		}
		f.closeRead();

		CurveFitter cf = new CurveFitter(x,y);
		System.out.println("starting doFit routine");
		cf.doFit(CurveFitter.IBEXFIT);
		}