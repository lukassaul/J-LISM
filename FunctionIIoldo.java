/**
* Utility class for passing around 2D functions
*
* Use this to get a 1D function by fixing one of the variables..
*  the index to fix is passed in ( 0 or 1 )
*/
public abstract class FunctionII {

	/**
	*  Override this for the function needed!!
	*
	*/
	public double function2(double x, double y) {
		return 0;
	}


	/**
	*
	* This returns a one dimensional function, given one of the values fixed
	*
	*/
	protected final Function getFunction(final int index, final double value) {
		if (index == 0) return new Function() {
			public double function(double x) {
				return function2(value,x);
			}
		};
		else if (index == 1) return new Function() {
			public double function(double x) {
				return function2(x,value);
			}
		};
		else {
			System.out.println("index out of range in FunctionII.getFunction");
			return null;
		}
	}
}
