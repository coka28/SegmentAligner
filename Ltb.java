package Alignment;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;

public class Ltb {
	
	String alignmentname;
	DataMaster DM;
	
	public Ltb(String name) {
		alignmentname = name;
		
	}
	
	public void loadDataMaster() throws ClassNotFoundException {
		String path = "C:/users/coka/desktop/programmierung/java/projekte/alignments/" + alignmentname + "/datamaster.object";
		
		try {
			FileInputStream fin = new FileInputStream(path);
			ObjectInputStream ois = new ObjectInputStream(fin);
			DM = (DataMaster)ois.readObject();
			ois.close();
		} catch (IOException e) {
			System.out.println("ERROR in loading DataMaster!");
			System.exit(0);
		}
	}
	
	public void getAll(double minscore) throws ClassNotFoundException, IOException {
		while (getOne(minscore));
		System.out.println("Found all alignments!");
	}
	
	public boolean getOne(double minscore) throws ClassNotFoundException, IOException {
		long globalindex = DM.getMaxIndex(minscore);
		if (globalindex==0) return false;
		Object[] tmp = DM.get(globalindex, true);
		double max = (double)tmp[0];
		boolean[] currentpointer = (boolean[])tmp[1];
		ArrayList<boolean[]> alignmentPointers = new ArrayList<boolean[]>();
		ArrayList<long[]> alignmentPositions = new ArrayList<long[]>();
		long[] endpos = DM.globalPosition(globalindex);
		
		ArrayList<Integer> touched = new ArrayList<Integer>();
		double newscore = 0;
		
		do {
			if (!touched.contains(DM.segmentIndex(globalindex))) touched.add(DM.segmentIndex(globalindex));
			DM.erase(globalindex);
			int zeros = 0;
			for (boolean i: currentpointer) if (!i) zeros++;
			alignmentPointers.add(currentpointer);
			alignmentPositions.add(DM.globalPosition(globalindex));
			long[] gpos = DM.globalPosition(globalindex);
			newscore += DM.match(gpos)-zeros*DM.pt;
			for (int i=0; currentpointer!=null && i<currentpointer.length; i++) if (currentpointer[i]) gpos[i]--;
			globalindex = DM.globalIndex(gpos);
			currentpointer = (boolean[]) DM.get(globalindex, true)[1];
		} while(currentpointer!=null);
		long[] startpos = DM.globalPosition(globalindex);
		for (int i=0; i<startpos.length; i++) startpos[i]++;
		
		String[] moddedstrings = new String[DM.strings.length];
		for (int i=0; i<moddedstrings.length; i++) moddedstrings[i]="";
		
		long[] currentpos = startpos;
		for (int i=alignmentPointers.size()-1; i>=0; i--) {
			currentpointer = alignmentPointers.get(i);
			currentpos = alignmentPositions.get(i);
			alignmentPointers.remove(i);
			alignmentPositions.remove(i);
			for (int k=0; k<currentpos.length; k++) {
				if (currentpointer[k]) {
					moddedstrings[k] += DM.strings[k].substring((int)currentpos[k], (int)currentpos[k]+1);
				} else moddedstrings[k] += "-";
			}
		}
		
		DM.reMax(touched);
		if (newscore<minscore) {
			return true;
		}
		print(moddedstrings, startpos, endpos, newscore);
		
		/**
		int cleanupstartfield = DM.segmentIndex(DM.globalIndex(startpos));
		int cleanupendfield = DM.segmentIndex(DM.globalIndex(endpos));
		DM.boxrescore(cleanupstartfield, cleanupendfield);*/
		//DM.linerescore(startpos, endpos);
		//DM.limited_rescore(present);
		
		if (newscore>14) DM.selective_rescore(touched);

		//DM.dumpAll(false);
		return true;
	}
	
	void print(String[] alignedStrings, long[] startpos, long[] endpos, double maxscore) throws IOException, ClassNotFoundException {
		BufferedWriter fw = new BufferedWriter(new FileWriter("C:/users/coka/desktop/programmierung/java/projekte/alignments/"+alignmentname+"/alignments.txt", true));
		
		for (int i=0; i<alignedStrings.length-1; i++)
			for (int k=i+1; k<alignedStrings.length; k++) {
				String[] tempstr = alignedStrings.clone();
				String str1 = "";
				String map = "";
				String str2 = "";
				int failures = 0;
				for (int n=0; n<tempstr[i].length(); n++) if (tempstr[i].charAt(n)==tempstr[k].charAt(n)) {
					if (tempstr[i].charAt(n)=='-') {
						tempstr[i] = tempstr[i].substring(0, n) + tempstr[i].substring(n+1);
						tempstr[k] = tempstr[k].substring(0, n) + tempstr[k].substring(n+1);
						n--;
					} else {
						str1 += tempstr[i].substring(n,n+1);
						map += "|";
						str2 += tempstr[k].substring(n,n+1);
					}
				} else {
					str1 += tempstr[i].substring(n,n+1);
					map += "*";
					str2 += tempstr[k].substring(n,n+1);
					failures++;
				}
				String newlines = "seq"+String.valueOf(i+1)+" "+String.valueOf(startpos[i])+".."+String.valueOf(endpos[i])+", "+
						"seq"+String.valueOf(k+1)+" "+String.valueOf(startpos[k])+".."+String.valueOf(endpos[k])+", score : "+
						String.valueOf(maxscore+"\n");
				int strlength = str1.length();
				int n=0;
				while (str1.length()>70*n) {
					String str1part = str1.substring(70*n, Math.min(strlength, 70*(n+1)));
					String str2part = str2.substring(70*n, Math.min(strlength, 70*(n+1)));
					String mappart = map.substring(70*n, Math.min(strlength, 70*(n+1)));
					newlines += str1part+"\n"+mappart+"\n"+str2part+"\n";
					n++;
				}
				newlines += "--Fehlerrate: "+String.valueOf((double)((10000*failures)/map.length())/100)+"%--\n";
				System.out.println(newlines);
				fw.write(newlines);
			}
			
			
		fw.close();}
	
}
