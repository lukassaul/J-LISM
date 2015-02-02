
public class TestFunction {
	public TestFunction() {
		MyF test = new MyF();
		MyF test1 = new MyF();
		MyF test2 = new MyF();
		MyF test3 = new MyF();

		Function a = test.getFunction(0,1.0);
		Function b = test1.getFunction(0,0.0);
		Function c = test2.getFunction(1,5.0);
		Function d = test3.getFunction(1,10.0);

		System.out.println(a.function(0.0)+" "+b.function(0.0)+" "+c.function(0.0)+" "+
				d.function(0.0));
	}

	public static final void main(String[] args) {
		TestFunction tf = new TestFunction();
	}
}


class MyF extends FunctionII {
	public double function2(double x, double y) {
		return Math.exp(x*x+y*y);
	}
}
