
//  THIS CLASS HAS SOME FILE UTILITY STUFF  //
//  Updated Version 02/28/00
//  Copyright 9 Elder
//	   nice synchronizations stuff here-
//new :  initRead, initWrite, closeRead, closeWrite, readLine, write
// improved:  readShitLn
////         be sure to use openFile() and closeFile() when using getLine and putLine
//			now catches those null pointers!
//
//  Does not support:  good error handling
//  				   politically correct java syntax

import java.util.*;
import java.io.*;
class file {
	private int numLines;
	private FileReader fr;
	private BufferedReader br;
	private FileWriter fw;
	private BufferedWriter bw;
	private boolean lineNumberCoherence;
	private boolean frozen = false;
	public String fileName = "";
	File test;  // used for functionality of Java's "File" class

	public file(String filename) {   // here's the constructor
		lineNumberCoherence = false;
		fileName = filename;
		numLines = 0;
	}

	public file() {
		lineNumberCoherence = false;
		fileName = "";
		numLines = 0;
	}

	private void o(String message) {
		System.out.println(message);
	}

	public synchronized void initRead() { // this is for manual file access
		o("initRead");
		frozen = true;
		try {
			fr = new FileReader(fileName);
			br = new BufferedReader(fr);
		}catch (Exception e) {
			System.out.println("Error- " + e.toString());
		}
	}

	public synchronized void initWrite(boolean append) {
		frozen = true;
		lineNumberCoherence=false;
		try {
			fw = new FileWriter(fileName, append);
			bw = new BufferedWriter(fw);
		}catch (Exception e) {
			System.out.println("Error- " + e.toString());
		}
	}

	public synchronized String readLine() {
		try {
			return br.readLine();
		}catch (Exception e) {
			System.out.println("Error- " + e.toString());
			return null;
		}
	}

	public synchronized void write(String shita) {
		try {
			bw.write(shita);
		}catch (Exception e) {
			System.out.println("Error- " + e.toString());
		}
	}

	public synchronized void closeRead() {
		frozen = false;
		try {
			br.close();
		}catch (Exception e) {
			System.out.println("Error- " + e.toString());
		}
	}

	public synchronized void closeWrite() {
		frozen = false;
		try {
			bw.close();
		}catch (Exception e) {
			System.out.println("Error- " + e.toString());
		}
	}


	// this should put shita (text) in a file
	public synchronized void saveShit (String shita, boolean append) {
		while (frozen) {};
		lineNumberCoherence = false;
		try {
			fw = new FileWriter(fileName,append);
			bw = new BufferedWriter(fw);
			bw.write(shita);
			bw.close();
		} catch (IOException e) {
			System.out.println("Error - " + e.toString());
		}
		catch (SecurityException b) {
			System.out.println("Go check security shit");
		}
		catch (NullPointerException c) {
			System.out.println("Null pointer, dude!"+c.toString());
		}
	}

	public synchronized void saveShit(String shita) {  // if you leave out the boolean, it's auto-overwrite
		saveShit(shita, false);
		}


	// 		this returns shita from a text file
	public synchronized String readShit() {
		checkFrozen();
		String textToBeReturned = null;
		try {
			String line = "";
			numLines = 0;
			fr = new FileReader(fileName);
			br = new BufferedReader(fr);
			boolean eof = false;
			while (!eof) {
			    line = br.readLine();
				if (line == null)
					eof = true;
				else {
					if (textToBeReturned == null) {
						textToBeReturned = line;
					}
					else {
						textToBeReturned += System.getProperty("line.separator")+line;
					}
					numLines++;
				}
			}
			br.close();
		} catch (IOException e) {
			System.out.println("Error - " + e.toString());
		}
		catch (SecurityException b) {
			System.out.println("Go check security shit");
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		return textToBeReturned;
	}

	// this is a little bit better- no tempFile means more memory taken up,
	// but we won't have any huge files here hopefully.
	public synchronized void deleteShitLn(int lineNumber) {
		checkFrozen();
		// we still have coherence...
		int currentLine=0;
		file oldFile = new file(fileName);
		String textToBeReturned = null;
		//file tempFile = new file("12345test.test");
		try {
			String line = "";
			fr = new FileReader(fileName);
			br = new BufferedReader(fr);
			boolean eof = false;
			while (!eof) {
				currentLine++;
			    line = br.readLine();
				if (line == null)
					eof = true;
				else if (currentLine != lineNumber) {
					if (textToBeReturned == null)
						textToBeReturned = line;
					else
						textToBeReturned += System.getProperty("line.separator")+line;
				}
			}
			br.close();
		}
		catch (IOException e) {
			System.out.println("Error - " + e.toString());
		}
		catch (SecurityException b) {
			System.out.println("Go check security shit");
		}
		catch (NullPointerException n) {
			System.out.println("Null Pointer in DeleteShitLn");
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		oldFile.saveShit(textToBeReturned);
		numLines--;
	}

	/*public synchronized void deleteShitLn(int lineNumber) {
		// we still have coherence...
		int currentLine=0;
		file oldFile = new file(fileName);
		String textToBeReturned = null;
		//file tempFile = new file("12345test.test");
		try {
			String line = "";
			char char1;
			char char2;
			FileReader letters = new FileReader(fileName);
			BufferedReader buff = new BufferedReader(letters);
			boolean eof = false;
			while (!eof) {
				currentLine++;
			    char1 = (char)buff.read();
			    char2 = (char)buff.read();
				if (char1 == null)
					eof = true;
				else if (currentLine != lineNumber) {
					if (textToBeReturned == null)
						textToBeReturned = line;
					else
						textToBeReturned += System.getProperty("line.separator")+line;
				}
			}
			buff.close();
		}
		catch (IOException e) {
			System.out.println("Error - " + e.toString());
		}
		catch (SecurityException b) {
			System.out.println("Go check security shit");
		}
		catch (NullPointerException n) {
			System.out.println("Null Pointer in DeleteShitLn");
		}
		oldFile.saveShit(textToBeReturned);
		numLines--;
	}*/

	public synchronized void changeShitLn(int lineNumber, String newStuff) {
		checkFrozen();
		int currentLine=0;
		file oldFile = new file(fileName);
		String textToBeReturned = null;
		//file tempFile = new file("12345test.test");
		try {
			String line = "";
			fr = new FileReader(fileName);
			br = new BufferedReader(fr);
			boolean eof = false;
			while (!eof) {
				currentLine++;
			    line = br.readLine();
				if (line == null)
					eof = true;
				else {
					if (currentLine != lineNumber) { //just add the line here
						if (textToBeReturned == null)
							textToBeReturned = line;
						else
							textToBeReturned += System.getProperty("line.separator")+line;
					}
					else {
						if (textToBeReturned == null) // here we add the new line
							textToBeReturned = newStuff;
						else
							textToBeReturned += System.getProperty("line.separator")+newStuff;
					}
				}
			}
			br.close();
		}
		catch (IOException e) {
			System.out.println("Error - " + e.toString());
		}
		catch (SecurityException b) {
			System.out.println("Go check security shit");
		}
		catch (NullPointerException n) {
			System.out.println("Null Pointer in ChangeShitLn");
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		oldFile.saveShit(textToBeReturned); // here we save it!
	}

	public synchronized String readShitLn(int lineNumber) {  // this returns the string from a line
		checkFrozen();
		String textToBeReturned = null;
		try {
			fr = new FileReader(fileName);
			br = new BufferedReader(fr);
			//LineNumberReader buff = new LineNumberReader(letters);
			String line = null;
			for (int i = 0; i<lineNumber; i++) {
				line = br.readLine();
			}
			//buff.setLineNumber(lineNumber+1);
			//textToBeReturned=buff.readLine();
			textToBeReturned = line;
			o("readShitLn got: " + textToBeReturned);
			br.close();
		} catch (IOException e) {
			System.out.println("Error - " + e.toString());
		} catch (SecurityException b) {
			System.out.println("Invalid security clearance- prepare to die");
		} catch (NullPointerException n) {
			System.out.println("Null pointer in readShitLn dude");
		}catch (Exception e) {
			e.printStackTrace();
		}
		return textToBeReturned;
	}


	// JON - these next routines go through the file 2N+1 times or some shit, compared to once above.
	// hopefully deprecated soon
	/* DO NOT DELETE*/
	public synchronized void jonDeleteShitLn(int lineNumber) {
		checkFrozen();
		file oldFile = new file(fileName);
		file tempFile = new file("12345test.test");
		int numLines = oldFile.readShitNumLines();
		tempFile.saveShit("");
		if (lineNumber <= numLines) {   // this method does nothing if linenumber > numlines
			for (int k = 1; k<lineNumber; k++) {
				tempFile.saveShit(oldFile.readShitLn(k)+System.getProperty("line.separator"),true);
				}
			if (numLines > lineNumber) {
				for (int j = lineNumber+1; j<numLines; j++) {
					tempFile.saveShit(oldFile.readShitLn(j)+
											System.getProperty("line.separator"),true);
					}
				tempFile.saveShit(oldFile.readShitLn(numLines),true); // no extra "\n"
				}
			oldFile.saveShit(tempFile.readShit() + "\n");
			boolean ahSo = tempFile.delete();
		}
		numLines--;
	}

	/* DO NOT DELETE */
	public synchronized void jonChangeShitLn(int lineNumber, String newStuff) {
		checkFrozen();
		file oldFile = new file(fileName);
		file tempFile = new file("12345test.test");
		int numLines = oldFile.readShitNumLines();
		tempFile.saveShit("");
		if (lineNumber <= numLines) {   // this method does nothing if linenumber > numlines
			for (int k = 1; k<lineNumber; k++) {
				tempFile.saveShit(oldFile.readShitLn(k)+System.getProperty("line.separator"),true);
				}
			if (lineNumber == numLines) {
				tempFile.saveShit(newStuff + "\n",true);
				}
			else {
				tempFile.saveShit(newStuff + System.getProperty("line.separator"), true);
				for (int j = lineNumber+1; j<numLines; j++) {
					tempFile.saveShit(oldFile.readShitLn(j)+
											System.getProperty("line.separator"),true);
					}
				tempFile.saveShit(oldFile.readShitLn(numLines),true);
			}
			oldFile.saveShit(tempFile.readShit() + "\n");
			tempFile.delete();
		}
	}

 	// use this for adding lines to a file in the fastest possible way
	public synchronized void addLines(int firstLineGoesHere, String[] lines) {
		// we still have line coherence, we know how many we are adding
		checkFrozen();
		System.out.println("starting file.addLInes-  first line:" + lines[0]);
		int numNewLines = lines.length;
		boolean tempIsEmpty = true;
		int currentLine=1; //this is the syntax.  First line = line 1 NOT line 0
		String line = "";
		//boolean done = false;
		try {
			// these are my boys
			fr = new FileReader(fileName);
			br = new BufferedReader(fr);
			String toBeReturned = "";

			boolean eof = false;
			while (!eof) {
				line=br.readLine();
				// first check if we need to add the goods and do it
				if(currentLine==firstLineGoesHere) {
					System.out.println("should be doing addLInes NOW");
					for(int i=0; i<numNewLines; i++) {
						if(tempIsEmpty) {
							toBeReturned = lines[i];
							tempIsEmpty = false;
						}
						else toBeReturned += System.getProperty("line.separator")+lines[i];
					}
					currentLine = currentLine + numNewLines;
				}
				// now check for eof
			    if (line == null)  eof = true;
			    else {
			    	if(tempIsEmpty) {
						toBeReturned = line;
						tempIsEmpty = false;
					}
					else toBeReturned += System.getProperty("line.separator")+line;
					currentLine++;
				}
			}
			saveShit(toBeReturned);
			numLines = numLines + numNewLines;
		}
		catch (IOException e) {
			System.out.println("Error - " + e.toString());
			}
		catch (SecurityException b) {
			System.out.println("Go check security shit");
		}
		catch (NullPointerException n) {
			System.out.println("Null Pointer in DeleteShitLn");
		}
		catch (Exception e) {
			System.out.println("Error - " + e.toString());
		}
	}

	public synchronized void addText(int firstLineGoesHere, String text) {
		checkFrozen();

		lineNumberCoherence = false;
		boolean tempIsEmpty = true;
		int currentLine=1; //this is the syntax.  First line = line 1 NOT line 0
		String toBeReturned = "";
		String line = "";
		boolean eof = false;
		try {
			// these are my boys
			fr = new FileReader(fileName);
			br = new BufferedReader(fr);
			while (!eof) {
				line=br.readLine();
				// first check if we need to add the goods and do it
				if(currentLine==firstLineGoesHere) {
					if(tempIsEmpty) {
						toBeReturned = text;
						tempIsEmpty = false;
					}
					else toBeReturned += System.getProperty("line.separator")+text;
				}
				// now check for eof
				if (line == null)  eof = true;
				else {
				   	if(tempIsEmpty) {
						toBeReturned = line;
						tempIsEmpty = false;
					}
					else toBeReturned += System.getProperty("line.separator")+line;
					currentLine++;
				}
			}
			saveShit(toBeReturned);
		}
		catch (IOException e) {
			System.out.println("Error - " + e.toString());
			}
		catch (SecurityException b) {
			System.out.println("Go check security shit");
		}
		catch (NullPointerException n) {
			System.out.println("Null Pointer in DeleteShitLn");
		}
		catch (Exception e) {
			System.out.println("Error - " + e.toString());
		}
	}


	public synchronized int readShitNumLines() {
		if (lineNumberCoherence) return numLines;
		String test = readShit(); //file could be changed or readShit not called yet.
		return numLines;
	}

	public synchronized boolean exists() {  // useful exists method
		test = new File(fileName);
		return test.exists();
	}

	public synchronized boolean delete() {  // useful delete method
		test = new File(fileName);
		return test.delete();
	}

	public synchronized long length() {  // returns size of file in bytes
		test = new File(fileName);
		return test.length();
	}
	public synchronized void setContents(file gehFile) { // this is quick!
		lineNumberCoherence = false;  // of a serious nature!
		File gf = new File(gehFile.fileName);
		test = new File(fileName);
		boolean niceAss = delete();
		try{
			gf.renameTo(test);
		}
		catch (Exception e) {
			System.out.println("Error in file class setCOntents: "+e.toString());
		}
	}
	private void checkFrozen() {
		if (frozen) {
			System.out.println("Important!  Somewhere there is no closeRead or closeWrite!!! "+fileName);
			try {Thread.sleep(1000);}
			catch(Exception e) {e.printStackTrace();}
		}
	}
}
// --- END FILE CLASS ---



