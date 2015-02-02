
/**
* Utility class for passing around 3D functions
*
* Use this to generate a 2D function by fixing one of the variables
*/
public abstract class FunctionIIII {

	/**
	*  Override this with an interesting 3D function
	*
	*/
	public double function4(double w, double x, double y, double z) {
		return 0;
	}


	/*
	* This returns a two dimensional function, given one of the values
	*
	*/
	public final FunctionIII getFunctionIII(final int index, final double value) {
		//value_ = value;
		if (index == 0) return new FunctionIII() {
			public double function3(double x, double y, double z) {
				return function4(value,x,y,z);
			}
		};
		else if (index == 1) return new FunctionIII() {
			public double function3(double x, double y, double z) {
				return function4(x,value,y,z);
			}
		};
		else if (index == 2) return new FunctionIII() {
			public double function3(double x, double y, double z) {
				return function4(x,y,value,z);
			}
		};
		else if (index == 3) return new FunctionIII() {
			public double function3(double x, double y, double z) {
				return function4(x,y,z,value);
			}
		};
		else {
			System.out.println("index out of range in FunctionIIII.getFunctionIII");
			return null;
		}
	}

	public static final void main(String[] args) {
	}
}

