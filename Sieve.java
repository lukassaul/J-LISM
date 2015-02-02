import java.lang.Math;
import java.util.Vector;
import java.util.Date;

public class Sieve {
	private static long MAX = (long)Math.pow(10,10);
	private static long MIN = (long)Math.pow(10,9);

	public long[] eNums;
	public int[] primes10000;
	private int numDigits = 20;
	private long numDivides = 0;
	public Sieve() {
		file ef = new file("e_digits.txt");
		ef.initRead();
		StringBuffer e = new StringBuffer("");
		eNums = new long[1000];
		int index = 0;
		String line="";
		while ((line=ef.readLine())!=null) {
			e.append(line);
		}
		System.out.println("read e file");
		for (int i=0; i<eNums.length; i++) {
			String sub = e.substring(i,i+numDigits);
			if (index < eNums.length-1)	eNums[index]=Long.parseLong(sub);
			System.out.println(eNums[index]+"");
			index++;
		}

		long num = MAX-MIN;

		System.out.println("num= " + num +" initializing vector.. ");
		Vector primes = new Vector();
		int debug = 0;
		for (int i=2; i<10000; i++) {
			primes.addElement(new Long(i));
		}
		System.out.println("finished init vector.. ");
		Date d1 = new Date();
		sieve(primes);
		Date d2 = new Date();
		System.out.println("sieve of 10000 took: " + (d2.getTime()-d1.getTime()));

		// sytem time test here..
		int doNothing;
		d1 = new Date();
		for (int i=0; i<1000000; i++) {
			doNothing = i/17;
		}
		d2 = new Date();
		System.out.println("10^6 divides took: " + (d2.getTime()-d1.getTime()));



		System.out.println("Found : " + primes.size() + "primes");
		primes10000 = new int[primes.size()];
		for (int i=0; i<primes.size(); i++) {
			primes10000[i] = ((Long)primes.elementAt(i)).intValue();
		}

		primes.clear();
		System.out.println("now we build our start vector");
		/*for (long l=MIN; l<MAX; l++) {
			debug ++;
			if (addIt(l,primes10000)) {
				primes.add(new Long(l));
			}
			if (debug==100000) {

				System.out.println("left: " + (num-l));

				System.out.println("primes.size: " + primes.size());

				debug = 0;
			}
		}*/

		// now we need to sieve..
		//	for (long l=101; l<
		d1 = new Date();
		long firstPrime = 0;
		int sequence = 0;
		int i = 0;
		while (i<eNums.length && firstPrime == 0) {
			if (checkPrime(eNums[i])) {
				if (firstPrime==0) {
					firstPrime=eNums[i];
					sequence = i;
				}
			}
			i++;
		}
		d2 = new Date();

		System.out.println(numDigits + " digits took: " + (d2.getTime()-d1.getTime()));
		System.out.println(sequence + " " + firstPrime);
		System.out.println("numDivides: " + numDivides);
	}

	private static boolean addIt(long l, long[] ll) {
		for (int i=0; i<ll.length; i++) {
			if (l%ll[i] == 0) return false;
		}
		return true;
	}

	/**
	*  The basic sieve -
	*   it starts from 2
	*   and goes to length of Vector..
	*
	*   not segmented!
	*/
	public static void sieve(Vector primes) {
		long tl=0; long ll=0;
		long min = ((Long)primes.elementAt(0)).longValue();
		if (min != 2) {
			System.out.println("incorrect use of sieve(v)");
			return;
		}
		for (int l=2; l<primes.size(); l++) {
			tl = ((Long)primes.elementAt(l)).longValue();
			int i=0;
			while (i<primes.size()) {
				ll=((Long)primes.elementAt(i)).longValue();
				if (tl>=ll) i++;
				else if (ll%tl==0) {
					// not a prime
					//System.out.println("not prime: " + ll);
					primes.remove(i);
				}
				else i++;
			}
		}
	}

	/**
	* works pretty well
	*
	*/
	public boolean checkPrime(long l) {
		// first check from our list of first few..
		for (int i=0; i<primes10000.length; i++) {
			numDivides++;
			if (l%primes10000[i]==0) {
				//System.out.println(l + " = " + primes10000[i] +" * " + l/primes10000[i]);
				return false;
			}
		}
		long s = (long)Math.sqrt((double)l);
		for (long i=3; i<=s; i+=2) {
			numDivides++;
			if (l%i==0) {
				System.out.println(l + " = " + i + " * " + l/i);
				return false;
			}
		}
		System.out.println(l + " is prime!");
		return true;
	}


	public static void main(String[] args) {
		Sieve a = new Sieve();

		// how many binaries < 10^10

	}
}