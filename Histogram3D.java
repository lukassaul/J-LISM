/**
*  Originally from CTOF dataset tools - Lukas Saul UNH '05
*
*/
public class Histogram3D {

	public int binsX, binsY, binsZ; // number of bins
	public int[][][] data; // this has a number of "events" per bin
	public float minX, minY, minZ, maxX, maxY, maxZ;
	public float widthX, widthY, widthZ;
	public float[] labelX; // this has the bottom value of each bin
	public float[] labelY;
	public float[] labelZ;
	public int totalEvents;

	public Histogram3D() {

		data = new int [0][0][0];
		labelX = new float[0];
		labelY = new float[0];
		labelZ = new float[0];
	}

	// a given histogram is made with one of these constructors
	public Histogram3D(float minX, float maxX, int _binsX,
						float minY, float maxY, int _binsY,
						float minZ, float maxZ, int _binsZ ) {
		//System.out.println("creating new histogram with numBins");
		binsX = _binsX;
		binsY = _binsY;
		binsZ = _binsZ;

		data = new int[binsX][binsY][binsZ];
		labelX = new float[binsX];
		labelY = new float[binsY];
		labelZ = new float[binsZ];

		widthX=(maxX-minX)/binsX;
		widthY=(maxY-minY)/binsY;
		widthZ=(maxZ-minZ)/binsZ;


		for (int i=0; i<binsX; i++) {
			labelX[i] = minX + i*(maxX-minX)/binsX;
		}
		for (int i=0; i<binsY; i++) {
			labelY[i] = minY + i*(maxY-minY)/binsY;
		}
		for (int i=0; i<binsZ; i++) {
			labelZ[i] = minZ + i*(maxZ-minZ)/binsZ;
		}

		for (int i=0; i<binsX; i++) {
			for (int j=0; j<binsY; j++) {
				for (int k=0; k<binsZ; k++) {
					data[i][j][k]=0;  // initiallize array to zero
				}
			}
		}
	}

	public void addEvent(float x, float y, float z) {

		if (minX>x | minY>y | minZ>z | maxX<x | maxY<y | maxZ < z) {
			System.out.println("Event not in histogram range");
			return;
		}

		int xBin = (int)((x-minX)/widthX);
		int yBin = (int)((y-minY)/widthY);
		int zBin = (int)((z-minZ)/widthZ);

		data[xBin][yBin][zBin]++;
		totalEvents++;
	}

	public Histogram2D getSlice(int axis, float min, float max)  {
		Histogram2D tbr;

		if (axis==0) {
			tbr = new Histogram2D(minY, maxY, binsY, minZ, maxZ, binsZ);

			// that also sets all the data to zeros, no events added

			if (min<minX | max>maxX) {
				System.out.println("wroXng coords for slice");
				return tbr;
			}
			int minBin = (int)((min-minX)/widthX);
			int maxBin = (int)((max-minX)/widthX);
			if (maxBin==minBin) {
				System.out.println("XSlice only 1 bin thick");
			}
			for (int j=0; j<binsY; j++) {
				for (int k=0; k<binsZ; k++) {
					for (int i=minBin; i<maxBin; i++) {
						tbr.data[j][k]+=data[i][j][k];  // just put all the counts in the bin directly
					}
				}
			}
		}
		else if (axis==1) {
			tbr = new Histogram2D(minX, maxX, binsX, minZ, maxZ, binsZ);
			if (min<minY | max>maxY) {
				System.out.println("wroYng coords for slice");
				return tbr;
			}
			int minBin = (int)((min-minY)/widthY);
			int maxBin = (int)((max-minY)/widthY);
			if (maxBin==minBin) {
				System.out.println("Y Slice only 1 bin thick");
			}
			for (int i=0; i<binsX; i++) {
				for (int k=0; k<binsZ; k++) {
					for (int j=minBin; j<maxBin; j++) {
						tbr.data[i][k]+=data[i][j][k];
					}
				}
			}
		}
		else if (axis==2) {
			tbr = new Histogram2D(minX, maxX, binsX, minY, maxY, binsY);
			if (min<minZ | max>maxZ) {
				System.out.println("wZrong coords for slice");
				return tbr;
			}
			int minBin = (int)((min-minZ)/widthZ);
			int maxBin = (int)((max-minZ)/widthZ);
			if (maxBin==minBin) {
				System.out.println("Z Slice only 1 bin thick");
			}
			for (int i=0; i<binsX; i++) {
				for (int j=0; j<binsY; j++) {
					for (int k=minBin; k<maxBin; k++) {
						tbr.data[i][j]+=data[i][j][k];
					}
				}
			}
		}
		else return tbr;
	}


}




