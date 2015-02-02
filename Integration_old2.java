import java.lang.Math;
import java.util.Date;

//Put your own function here and intgrate it
// Lukas Saul - Warsaw 2000

// lets use someone elses, shall we?
// nah, they all suck.  I'll write my own


public class Integration {
	private int startDivision = 5;
	private long maxDivision = 100000;
	private int maxTime = 60;
	private boolean checkTime = false;
	private double maxError = .1;
	private double averageValue;
	private Function f;
	//private int splitSize;
	private double f1, f2, f3, f4, f5, f6, f7;
	private Date d1, d2, d3;
	private double ll,ul,last;

	public void setStartDivision(int s) { startDivision = s; }
	public int getStartDivision() {return startDivision; }
	public void setMaxDivision(long m) {maxDivision = m; }
	public long getMaxDivision() {return maxDivision;}
	public void setMaxTime(int mt) {maxTime = mt;}
	public int getMaxTime() {return maxTime; }
	public void setCheckTime(boolean ct) {checkTime = ct; }
	public boolean getCheckTime() {return checkTime; }
	public void setMaxError(double me) {	if (me<1 & me>0) maxError = me; }
	public double getMaxError() {return maxError;}
	public void setFunction(Function fi) {f=fi;}
	//public void setSplitSize (int ss) {splitSize = ss;}
	//public int getSplitSize () {return splitSize; }


	private double IntegrateSimple(double lowerLimit, double upperLimit) {
		if (lowerLimit==ll & upperLimit == ul) return last;
		else {
			ll = lowerLimit; ul=upperLimit;
			f1 = f.function(lowerLimit); f2 = f.function(upperLimit);
			last = (upperLimit - lowerLimit) * (f1 + .5*(f2-f1)); // trapezoidal area
			return last;
		}
	}

	public double integrate (double lowerLimit, double upperLimit) {
		if (checkTime) d1 = new Date();
		f3 = 0;
		f4 = (upperLimit - lowerLimit)/startDivision; // f4 = divisionSize
		for (int i = 0; i< startDivision; i++) {
			f3 += recurse(lowerLimit+(i*f4), lowerLimit+((i+1)*f4), startDivision);
		}
		return f3;
	}

	private double recurse (double lower, double upper, long depth) {
		if (checkTime) {
			d2 = new Date();
			if ((d2.getTime() - d1.getTime())/1000 > maxTime) {
				return IntegrateSimple(lower, upper);
			}
		}

		if (depth >= maxDivision) return IntegrateSimple(lower, upper);

		f5 = IntegrateSimple(lower, lower + (upper-lower)/2);
		f5 += IntegrateSimple(lower+ (upper-lower)/2, upper);
		f6 = IntegrateSimple(lower, upper);

		// if the difference between these two intgrals is within error, return the better
		if (Math.abs(1-f6/f5) < maxError) return f5;
		// if not divide up the two and recurse
		else return recurse(lower, lower + (upper-lower)/2, depth*2) +
					recurse(lower+ (upper-lower)/2, upper, depth*2);

	}
}
