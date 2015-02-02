import java.util.StringTokenizer;
import java.util.Random;
public class Benford2 {
	public file inFile, outFile;
	public int[][] digitTallies;
	//public double[] factors = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0}; // last one placeholder
	public int skipLines, dataColumn;
	public StringTokenizer st;
	public Random r;

	public Benford2(String[] args) {
		r=new Random();
		try {
			skipLines = Integer.parseInt(args[1]);
			dataColumn = Integer.parseInt(args[2]);
			digitTallies = new int[10][2];

			o("datacolumn: " + dataColumn + " skiplines: " + skipLines );

			// set all digit tallies to zero
			for (int i=0; i<2; i++) {
				for (int j=0; j<10; j++) {
					digitTallies[j][i]=0;
				}
			}

			inFile = new file(args[0]);
			inFile.initRead();
			String line = "";
			for (int i=0; i<skipLines; i++) line=inFile.readLine();
			o("made it here 2");

			double max = -Double.MAX_VALUE;
			double min = Double.MAX_VALUE;
			// lets go through and find max and min of the data
			double avg = 0.0;
			int index = 0;
			while ((line=inFile.readLine())!=null) {
				String numString = "";
				st = new StringTokenizer(line);
				for (int j=0; j<dataColumn; j++) {
					numString = st.nextToken();
				}
				double theNumber = Double.parseDouble(numString);
				if (theNumber!=-1) {
					if (theNumber<min) min=theNumber;
					if (theNumber>max) max=theNumber;
					avg+=theNumber;
					index++;
				}
			}
			avg/=(double)index;
			inFile.closeRead();
			o("max: " + max);
			o("min: " + min);
			o("avg: " + avg);

			inFile = new file(args[0]);
			inFile.initRead();
			for (int i=0; i<skipLines; i++) line=inFile.readLine();
			while ((line=inFile.readLine())!=null) {
				String numString = "";
				st = new StringTokenizer(line);
				for (int j=0; j<dataColumn; j++) {
					numString = st.nextToken();
				}
				double theNumber = Double.parseDouble(numString);
				if (theNumber != -1) {

					//o("numstring: " + numString);
					for (int i=0; i<2; i++) {
						if (i==1) theNumber = (theNumber-avg);
						if (theNumber!=0.0) {
							//o("thenumber " + theNumber);
							int theDigit = getDigit(theNumber);
							digitTallies[theDigit][i]++;
						}
					}

				}
			}
			o("made it here 3");
			// we have the tallies -- lets generate a nice data file
			outFile = new file("benford_results2.txt");
			outFile.initWrite(false);
			for (int j=0; j<10; j++) {
				line = "";
				for (int i=0; i<2; i++) {
					line +=	digitTallies[j][i]+"\t";
				}
				outFile.write(line+"\n");
			}
			outFile.closeWrite();
			// done?
		}
		catch (Exception e) {
			o("Format:  java Benford filename.ext numLinesToSkip dataColumn_start_with_1");
			e.printStackTrace();
		}
	}

	public int getDigit(double d) {
		if (d<0) d*=-1.0;
		while (d<=1.0) d*=10.0;
		while (d>=10.0) d/=10.0;
		return (int)(Math.floor(d));
	}

	public static final void main(String[] args) {
		Benford2 b = new Benford2(args);
	}

	public static void o(String s) {
		System.out.println(s);
	}
}