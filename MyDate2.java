

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.text.DefaultEditorKit;
import java.net.*;
import java.io.*;
import java.util.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

// Now keeps the last dates in memory to avoid parsing again.
// Also remembers most recent format and uses it, rather than trying with each one.
public class MyDate2 {
	private int lastPattern=0;
	private String lastDateString="";
	private Date lastDate=null;
	private Date toBeReturned=null;
	private boolean didIt = false;
	private int iCounter = 0;
	private SimpleDateFormat sdf;
	private String pattern[] = {
					"EEE, dd MMM yyyy HH:mm:ss z", // this one works for most  services
					"EEE, dd MMM yyyy HH:mm:ss zzzz", // does this help?
					"dd MMM yyyy HH:mm:ss z", // this one works for hotmail
					"EEE MMM dd HH:mm:ss zzz yyyy",
					"d MMM yy HH:mm:ss zzz",
					"EEE, dd MMM yyyy HH:mm:ss", // Needed to aid at hotmail
					"dd MMM yyyy HH:mm:ss", // Another hotmail fix: 7 Feb 2000 01:30:18
					"EEE, dd MMM yy hh:mma z", //usa.net?
					"EEE, dd MMM yy HH:mm:ss z", // Hotmail strikes: Sat, 12 Aug 00 12:03:18 Pacific Daylight Time
					"dd MMM yy hh:mm:ss a", // Hotmail strikes: 11 Aug 00 7:29:09 AM
					"EEE, dd MMM yy HH:mma z", // Usa.net: Thu, 21 Sep 00 13:18PM PDT
					"EEEE, MMMM d, yyyy h:mm:ss 'o''clock' a z", // the rest here are guesses
					"MMMM d, yyyy hh:mm:ss a z",
					"dd-MMM-yy h:mm:ss a",
					"MM/dd/yy h:mm a",
					"EEE, dd MMM yyyy HH:mm z", /// Tue, 19 Dec 2000 15:57 +0000
					"EEE, dd MMM yyyy :mm z", /// Thu, 2 Nov 2000 :21 EST
					"dd-MM-yy hh:mm:ss a",
					"MM/dd/yy HH:mm:ss",
					"EEE, MMM dd, yyyy hh:mm:ss a", // the new netscape: Saturday, January 13, 2001 1:47:58 AM
					"EEE, dd MMM yyyy HH:mm", // Australia, telstra.com: Sunday, 7 January 2001 23:16
					"HH:mm:ss MM/dd/yy" };
	//NOTE KEY HERE: under 4 letters - look for abbreviation. 4 0r more for written out.
	// H - 24 hr scale h - 12 hr scale a - am/pm E - day of week - see API
/*

 Symbol   Meaning                 Presentation        Example
 ------   -------                 ------------        -------
 G        era designator          (Text)              AD
 y        year                    (Number)            1996
 M        month in year           (Text & Number)     July & 07
 d        day in month            (Number)            10
 h        hour in am/pm (1~12)    (Number)            12
 H        hour in day (0~23)      (Number)            0
 m        minute in hour          (Number)            30
 s        second in minute        (Number)            55
 S        millisecond             (Number)            978
 E        day in week             (Text)              Tuesday
 D        day in year             (Number)            189
 F        day of week in month    (Number)            2 (2nd Wed in July)
 w        week in year            (Number)            27
 W        week in month           (Number)            2
 a        am/pm marker            (Text)              PM
 k        hour in day (1~24)      (Number)            24
 K        hour in am/pm (0~11)    (Number)            0
 z        time zone               (Text)              Pacific Standard Time
 '        escape for text         (Delimiter)
 ''       single quote            (Literal)           '

The count of pattern letters determine the format.
(Text): 4 or more pattern letters--use full form, < 4--use short or abbreviated form if one exists.
(Number): the minimum number of digits. Shorter numbers are zero-padded to this amount.
Year is handled specially; that is, if the count of 'y' is 2, the Year will be truncated to 2 digits.
(Text & Number): 3 or over, use text, otherwise use number.
*/
	//here's the main routine:
	public Date parse(String date) {

		// In order to be a little efficient, we are going to try
		// and check our last known good date against this one
		// Some services standardize on dates, so we'll be able
		// to be more efficient and then faster

		if (date.equals(lastDateString)) return lastDate;
		// if so, we audi.

		sdf = new SimpleDateFormat(pattern[lastPattern]);
		sdf.setLenient(false); // lenient my ass!
		try {
			toBeReturned = sdf.parse(date);
			didIt = true;
			lastDateString = date;
			lastDate = toBeReturned;
		}
		catch (Exception e) {didIt = false; }
		if (didIt) return toBeReturned;
		// if so, we audi again

		// counter needs to be reset so that we start from the beginning of our patterns
		iCounter = 0;

		// worst case scenario, loop through all patterns
		while (!didIt) {
		//	System.out.println("Checking pattern: " + pattern[iCounter]);
			sdf.applyPattern(pattern[iCounter]);
			sdf.setLenient(false); // lenient my ass!

		//	System.out.println("Checking pattern: " + pattern[iCounter]);

			try {
				toBeReturned = sdf.parse(date);
				didIt = true;
				lastDateString = date;
				lastDate = toBeReturned;
			}
			catch (Exception e) {didIt = false; }
			if (didIt) return toBeReturned;
			// we audi
	//		System.out.println("date: " + date);

			if (iCounter==pattern.length-1) { // so we couldn't parse it!
				System.out.println("FIX THIS!!  UNABLE TO PARSE: " + date);
				toBeReturned = new Date();
				return toBeReturned;
			}
			else iCounter++;
		}
		// We should never get here
		System.out.println("problems in the date parser of a serious nature!!");
		return toBeReturned; //required for compilation
	}
}
