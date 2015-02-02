import java.util.StringTokenizer;
import java.util.Date;
import java.util.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;


/**
*  Ephemeris info here!
*  Position of Earth, pointing of IBEX, position of IBEX based on DOY input
*
*  Orbit info also
*/
public class EarthIBEX {

	// read these in J2000 coordinates
	public double[] rx,ry,rz;
	public double[] vx,vy,vz;
	public Date[] orbStarts, orbEnds;
	public Date[] dates;
	public HelioVector[] pointings;
	public Date startDate;
	public long currentDate;
	public file inFile;
	public file orbFile = new file("orbit_times.txt");
	public StringTokenizer st;
	public int nn = 24*120; // total amount of entries in earth info file, 1 per hour for 120 days
	public int norbs = 120;  // number of orbits in orbit times file
	public TimeZone tz = TimeZone.getTimeZone("UTC");
	public Calendar c;
	public SimpleDateFormat sdf;
	public String line;
	public EarthIBEX(int year) {

		// read the orbit times
		orbStarts = new Date[norbs];
		orbEnds = new Date[norbs];
		orbFile.initRead();
		String pattern = "yyyy/MM:dd:HH:ss z";
		sdf = new SimpleDateFormat(pattern);
		sdf.setLenient(false);
		//String test = "2010/06:26:17:44";
		//try {
		//	Date dt = sdf.parse(test + " UTC");
		//	System.out.println("should be 2010/06:26:17:44 " + dt.toGMTString());
		//}
		//catch (Exception e) {
		//	e.printStackTrace();
		//}
		try {
			line = orbFile.readLine();  // header
			for (int i=0; i<norbs; i++) {
				line = orbFile.readLine();
				st = new StringTokenizer(line);
				String garbage = st.nextToken();
				orbStarts[i] = sdf.parse(st.nextToken()+ " UTC");
				System.out.println("read orbit date: " + i + " = " +orbStarts[i].toGMTString());
			}
			orbFile.closeRead();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		// done reading orbit times



		// now lets read the pointing information from the ECLIPJ2000 idl output
		pointings = new HelioVector[100];
		System.out.println(" trying to read pointing information ");
		file f = new file("ibex_point.txt");
		f.initRead();
		line = f.readLine();
		line = f.readLine(); // throw away header
		for (int i=6; i<93; i++) {
			line = f.readLine();
			st = new StringTokenizer(line);
			int oo = (int)Double.parseDouble(st.nextToken());
			double phi = Double.parseDouble(st.nextToken());
			double theta = Double.parseDouble(st.nextToken());
			System.out.println("pointing for " + oo + " " + phi + " " + theta);
			pointings[oo]=new HelioVector(HelioVector.SPHERICAL,1.0,phi*Math.PI/180.0,theta*Math.PI/180.0);
		}
		f.closeRead();
		// done reading pointing info!!


		// read the earth pos. and velocity data from the file made by IDL
		rx = new double[nn];
		ry = new double[nn];
		rz = new double[nn];
		vx = new double[nn];
		vy = new double[nn];
		vz = new double[nn];
		dates = new Date[nn];
		if (year==2009) {
			inFile = new file("earth_pos_09.txt");
			c = Calendar.getInstance();
			c.setTimeZone(tz);
			c.set(Calendar.YEAR, 2009);
			c.set(Calendar.MONTH, 1);
			c.set(Calendar.DAY_OF_YEAR, 1);
			c.set(Calendar.HOUR_OF_DAY, 0);
			c.set(Calendar.MINUTE, 0);
			c.set(Calendar.SECOND, 0);
			startDate = c.getTime();
			System.out.println("start date: " +  startDate.toString());
			currentDate = startDate.getTime();
		}
		if (year==2010) {
			inFile = new file("earth_pos_10.txt");
			c = Calendar.getInstance();
			c.setTimeZone(tz);
			c.set(Calendar.YEAR, 2010);
			c.set(Calendar.MONTH, 1);
			c.set(Calendar.DAY_OF_YEAR, 1);
			c.set(Calendar.HOUR_OF_DAY, 0);
			c.set(Calendar.MINUTE, 0);
			c.set(Calendar.SECOND, 0);
			startDate = c.getTime();
			System.out.println("start date: " +  startDate.toString());
			currentDate = startDate.getTime();
		}
		inFile.initRead();
		line = "";
		// skip the header
		for (int i=0; i<2; i++) line=inFile.readLine();
		// read all the data

		// these are hourly datas
		for (int i=0; i<nn; i++) {
			dates[i]=new Date(currentDate);
			line = inFile.readLine(); // garbage
			line = inFile.readLine(); // r
			st = new StringTokenizer(line);
			rx[i]=Double.parseDouble(st.nextToken());
			ry[i]=Double.parseDouble(st.nextToken());
			rz[i]=Double.parseDouble(st.nextToken());
			line = inFile.readLine(); // v
			st = new StringTokenizer(line);
			vx[i]=Double.parseDouble(st.nextToken());
			vy[i]=Double.parseDouble(st.nextToken());
			vz[i]=Double.parseDouble(st.nextToken());
			currentDate+=3600*1000;
		}
	}

	/**
	* Return the Earth's position on that day in Heliospheric coords
	*
	*  could be much faster with a hashtable but this is not our limiting routine
	*  so we do it the hack way
	*/
	public HelioVector getEarthPosition(Date d) {
		long dd = d.getTime();
		int ourIndex = 0;
		for (int i=0; i<nn; i++) {
			if (dd<dates[i].getTime()) {
				ourIndex=i;
				i=nn;
			}
		}
		return new HelioVector(HelioVector.CARTESIAN, rx[ourIndex]*1000, ry[ourIndex]*1000, rz[ourIndex]*1000);
	}


	/**
	* Return the point direction of IBEX spin axis.  This should be toward the sun
	*
	*/
	public HelioVector getIbexPointing(Date d) {
		int orbit1 = getOrbitNumber(d);
		return getIbexPointing(orbit1);
	}

	/**
	* Return the point direction of IBEX spin axis.  This should be toward the sun
	*
	*/
	public HelioVector getIbexPointing(int orb) {
		if (orb<6 | orb>93) {
			System.out.println(" IbexPointing not available!! ");
			return new HelioVector();
		}
		return pointings[orb];
	}


	/**
	*  Return the velocity vector of Earth moving about the sun
	*
	*/
	public HelioVector getEarthVelocity(Date d) {
		long dd = d.getTime();
		int ourIndex = 0;
		for (int i=0; i<nn; i++) {
			if (dd<dates[i].getTime()) {
				ourIndex=i;
				i=nn;
			}
		}
		return new HelioVector(HelioVector.CARTESIAN, vx[ourIndex]*1000, vy[ourIndex]*1000, vz[ourIndex]*1000);
	}

	/**
	* Get the orbit number for a date
	*
	*/
	public int getOrbitNumber(Date d) {
		long dd = d.getTime();
		int ourIndex = 0;
		for (int i=0; i<nn-1; i++) {
			if (dd>orbStarts[i].getTime() & dd<orbStarts[i+1].getTime()) {
				ourIndex=i+1;
				i=nn;
			}
		}
		//System.out.println(" orbit number for " +d.toString()+" is " + ourIndex);
		return ourIndex;
	}


	public Date getDate(double doy, int year) {
		int day = (int)doy;
		double hour = (doy-day)*24.0;
		//System.out.println("doy: " + doy + " day+ " + day + " hour " + hour);
		c = Calendar.getInstance();
		c.setTimeZone(tz);
		c.set(Calendar.YEAR, year);
		c.set(Calendar.DAY_OF_YEAR, (int)doy);
		c.set(Calendar.HOUR_OF_DAY, (int)(hour));
		c.set(Calendar.MINUTE, getMinute(hour));
		c.set(Calendar.SECOND, getSecond(hour));
		//352.97629630  2008 = 913591566
		Date d = c.getTime();
		return d;
	}

	public static int getMinute(double hour) {
		double minutes = hour*60.0;
		// subtract hours
		double iHour = (double)(int)(hour)*60.0;
		minutes = minutes - iHour;
		//System.out.println("hour " + hour + " minutes" + minutes);
		return (int)(minutes);
	}

	public static int getSecond(double hour) {
		double seconds = hour*3600.0;
		double iHour = (int)(hour)*3600.0;
		double iMinute = getMinute(hour)*60.0;
		// subtract hours and minutes
		seconds = seconds - iHour - iMinute;
		//System.out.println("rem: "+ rem);
		return (int)(seconds);
	}

	/**
	* For testing
	*/
	public static final void main(String[] args) {
		try {
			String pattern = "yyyy/MM:dd:HH:ss z";
			SimpleDateFormat sdf2 = new SimpleDateFormat(pattern);
			sdf2.setLenient(false);

			EarthIBEX ei = new EarthIBEX(2009);
			Date d = sdf2.parse("2009/01:31:23:44 UTC");
			System.out.println("earth pos: " + ei.getEarthPosition(d).toAuString());
			System.out.println("earth vel: " + ei.getEarthVelocity(d).toKmString());
			System.out.println("ibex point: " + ei.getIbexPointing(d).toString());

			// now we output position and pointing (longitude j2000 ecliptic)  (day 10 - 120)
			file ff = new file("fig_1_2009.txt");
			ff.initWrite(false);
			double step = 1.0/24.0;  // step by 1 hour
			for (double doy=10.0; doy<120; doy+=step) {
				HelioVector pos = ei.getEarthPosition(ei.getDate(doy,2009));
				HelioVector point = ei.getIbexPointing(ei.getDate(doy,2009));
				ff.write(doy + "\t" + pos.getPhi()*180.0/Math.PI + "\t" + point.getPhi()*180.0/Math.PI+"\n");
			}
			ff.closeWrite();



			// test get date
			System.out.println(ei.getDate(2.1, 2009).toString());
			System.out.println(ei.getDate(2.15, 2010).toString());
			System.out.println(ei.getDate(3.5, 2010).toString());
			System.out.println(ei.getDate(2.6, 2010).toString());
			System.out.println(ei.getDate(2.7, 2010).toString());
			System.out.println(ei.getDate(2.8, 2010).toString());

		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}


}


/*
We need to keep track of Earth's Vernal Equinox
and use J2000 Ecliptic Coordinates!!!

March
2009 	20 	11:44
2010 	20 	17:32
2011 	20 	23:21
2012 	20 	05:14
2013 	20 	11:02
2014 	20 	16:57
2015 	20 	22:45
2016 	20 	04:30
2017 	20 	10:28

This gives the location of the current epoch zero ecliptic longitude..
However


/*

/*
Earth peri and aphelion

		perihelion			aphelion
2007 	January 3 	20:00 	July 7 	00:00
2008 	January 3 	00:00 	July 4 	08:00
2009 	January 4 	15:00 	July 4 	02:00
2010 	January 3 	00:00 	July 6 	12:00
2011 	January 3 	19:00 	July 4 	15:00
2012 	January 5 	01:00 	July 5 	04:00
2013 	January 2 	05:00 	July 5 	15:00
2014 	January 4 	12:00 	July 4 	00:00
2015 	January 4 	07:00 	July 6 	20:00
2016 	January 2 	23:00 	July 4 	16:00
2017 	January 4 	14:00 	July 3 	20:00
2018 	January 3 	06:00 	July 6 	17:00
2019 	January 3 	05:00 	July 4 	22:00
2020 	January 5 	08:00 	July 4 	12:00

*/