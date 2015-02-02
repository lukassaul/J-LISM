/**
*  Originally from CTOF dataset tools - Lukas Saul UNH '05
*
*/
public class Histogram2D {

	public int binsX, binsY; // number of bins
	public int[][] data; // this has a number of "events" per bin
	public float minX, minY, maxY, maxX;
	public float widthX, widthY;
	public float[] labelX; // this has the bottom value of each bin
	public float[] labelY;
	public int totalEvents;

	public Histogram2D() {

		data = new int [0][0];
		labelX = new float[0];
		labelY = new float[0];
	}

	// a given histogram is made with one of these constructors
	public Histogram2D(float minX, float maxX, int _binsX,
						float minY, float maxY, int _binsY) {
		//System.out.println("creating new histogram with numBins");
		binsX = _binsX;
		binsY = _binsY;

		data = new int[binsX][binsY];
		labelX = new float[binsX];
		labelY = new float[binsY];

		widthX=(maxX-minX)/binsX;
		widthY=(maxY-minY)/binsY;


		for (int i=0; i<binsX; i++) {
			labelX[i] = minX + i*(maxX-minX)/binsX;
		}
		for (int i=0; i<binsY; i++) {
			labelY[i] = minY + i*(maxY-minY)/binsY;
		}

		for (int i=0; i<binsX; i++) {
			for (int j=0; j<binsY; j++) {
					data[i][j]=0;  // initiallize array to zero
			}
		}
	}

	public void addEvent(float x, float y) {

		if (minX>x | minY>y  | maxX<x | maxY<y ) {
			System.out.println("Event not in histogram range");
			return;
		}

		int xBin = (int)((x-minX)/widthX);
		int yBin = (int)((y-minY)/widthY);

		data[xBin][yBin]++;
		totalEvents++;
	}

	public void setSlice(int axis, float min, float max)  {

		if (axis==0) {

		}
	}

	public double[] getArray() {
		double[] tbr = new double[data.length*data[0].length];
		int counter = 0;
		for (int i=0; i<data.length; i++) {
			for (int j=0; j<data[i].length; j++) {
				tbr[counter] = (double)data[i].data[j];
				counter++;
			}
		}
		return tbr;
	}



}




