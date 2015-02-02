import java.util.StringTokenizer;
import java.util.Random;
public class Benford {
	public file inFile, outFile;
	public int[][] digitTallies;
	public double[] factors = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0}; // last one placeholder
	public int skipLines, dataColumn;
	public StringTokenizer st;
	public Random r;

	public Benford(String[] args) {
		r=new Random();
		try {
			skipLines = Integer.parseInt(args[1]);
			dataColumn = Integer.parseInt(args[2]);
			digitTallies = new int[10][factors.length];

			o("datacolumn: " + dataColumn + " skiplines: " + skipLines );

			// set all digit tallies to zero
			for (int i=0; i<factors.length; i++) {
				for (int j=0; j<10; j++) {
					digitTallies[j][i]=0;
				}
			}

			inFile = new file(args[0]);
			inFile.initRead();
			String line = "";
			for (int i=0; i<skipLines; i++) line=inFile.readLine();
			o("made it here 2");

			while ((line=inFile.readLine())!=null) {
				String numString = "";
				st = new StringTokenizer(line);
				for (int j=0; j<dataColumn; j++) {
					numString = st.nextToken();
				}
				//o("numstring: " + numString);
				for (int i=0; i<factors.length; i++) {
					double theNumber = Double.parseDouble(numString);
					if (theNumber!=0.0) {
						//o("thenumber " + theNumber);
						if (i==factors.length-1) factors[9]=1+r.nextInt(9);
						theNumber = theNumber*factors[i];
						int theDigit = getDigit(theNumber);
						digitTallies[theDigit][i]++;
					}
				}
			}
			o("made it here 3");
			// we have the tallies -- lets generate a nice data file
			outFile = new file("benford_results.txt");
			outFile.initWrite(false);
			for (int j=0; j<10; j++) {
				line = "";
				for (int i=0; i<factors.length; i++) {
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
		while (d<=1.0) d*=10.0;
		while (d>=10.0) d/=10.0;
		return (int)(Math.floor(d));
	}

	public static final void main(String[] args) {
		Benford b = new Benford(args);
	}

	public static void o(String s) {
		System.out.println(s);
	}
}