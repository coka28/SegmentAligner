package Alignment;

import java.io.BufferedReader;  
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class GlobAl {
	
	Scoring S;
	
	public GlobAl(String[] strings, float match, float sub, float penalty) throws IOException {
		S = new Scoring(strings);
		S.mt = match;
		S.mm = sub;
		S.pt = penalty;
		S.score(new boolean[0]);
		tb();
	}
	
	public void tb() throws IOException {

	int index = S.factorials[S.strings.length-1]*S.strLengths[S.strings.length-1]-1;
	
	boolean[] currentpointer = S.pointers[index];
	ArrayList<boolean[]> alignmentPointers = new ArrayList<boolean[]>();
	ArrayList<int[]> alignmentPositions = new ArrayList<int[]>();
	while(index != 0) {
		alignmentPointers.add(currentpointer);
		alignmentPositions.add(S.position(index));
		for (int i=0; i<currentpointer.length; i++) if (currentpointer[i]) index -= S.factorials[i];
		currentpointer = S.pointers[index];
	}
		
	String[] moddedstrings = new String[S.strings.length];
	for (int i=0; i<moddedstrings.length; i++) moddedstrings[i]="";
	
	int[] currentpos = S.position(S.factorials[S.strings.length-1]*S.strLengths[S.strings.length-1]-1);
	System.out.println(alignmentPointers.size());
	for (int i=alignmentPointers.size()-1; i>=0; i--) {
		currentpointer = alignmentPointers.get(i);
		currentpos = alignmentPositions.get(i);
		alignmentPointers.remove(i);
		alignmentPositions.remove(i);
		for (int k=0; k<currentpos.length; k++) {
			if (currentpointer[k]) {
				moddedstrings[k] += S.strings[k].substring((int)currentpos[k], (int)currentpos[k]+1);
			} else moddedstrings[k] += "-";
		}
	}
	
	
	long[] startpos = new long[S.strings.length];
	long[] endpos = new long[S.strings.length];
	print(moddedstrings, startpos, endpos);
	
	}

	void print(String[] alignedStrings, long[] startpos, long[] endpos) throws IOException {
		BufferedWriter fw = new BufferedWriter(new FileWriter("globalalignment.txt"));
		
		for (int i=0; i<alignedStrings.length-1; i++)
			for (int k=i+1; k<alignedStrings.length; k++) {
				String[] tempstr = alignedStrings;
				String str1 = "";
				String map = "";
				String str2 = "";
				for (int n=0,m=0; n<tempstr[i].length() && m<tempstr[k].length(); m++,n++) if (tempstr[i].charAt(n)==
						tempstr[k].charAt(m)) {
					if (tempstr[i].charAt(n)=='-') {
						tempstr[i] = tempstr[i].substring(0, n) + tempstr[i].substring(n+1);
						tempstr[k] = tempstr[k].substring(0, n) + tempstr[k].substring(n+1);
						n--;
						m--;
					} else {
						str1 += tempstr[i].substring(n, n+1);
						map += "|";
						str2 += tempstr[k].substring(n, n+1);
					}
				} else {
					str1 += tempstr[i].substring(n, n+1);
					map += "*";
					str2 += tempstr[k].substring(n, n+1);
				}
				String newlines = "str"+String.valueOf(i+1)+" : ["+String.valueOf(startpos[i])+", "+String.valueOf(endpos[i])+"[ ; "+
						"str"+String.valueOf(k+1)+" : ["+String.valueOf(startpos[k])+", "+String.valueOf(endpos[k])+"[\n";
				int strlength = str1.length();
				int n=0;/**
				while (str1.length()>0) {
					String str1part = str1.substring(70*n, Math.min(strlength-1, 70*(n+1)));
					str1 = str1.substring(Math.min(strlength-1, 70*(n+1)));
					String str2part = str2.substring(70*n, Math.min(strlength-1, 70*(n+1)));
					str2 = str2.substring(Math.min(strlength-1, 70*(n+1)));
					String mappart = map.substring(70*n, Math.min(strlength-1, 70*(n+1)));
					map = map.substring(Math.min(strlength-1, 70*(n+1)));
					newlines += str1part+"\n"+mappart+"\n"+str2part+"\n";
					n++;
				}*/
				newlines += str1+"\n"+map+"\n"+str2+"\n";
				newlines += "\n";
				System.out.println(newlines);
				fw.write(newlines);
			}
			
			
		fw.close();}
}