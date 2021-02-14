package Alignment;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;

public class LScore {
	
	public String[] strings;
	int[] partdims;
	public int[] partsize;
	int[] factorials;
	long totalpartsize = 1;
	int totalpartnum = 1;
	double[] maxmatrix;
	int[] maxindexmatrix;
	String alignmentname;
	double mt=1, mm=-9, pt=9;
	
	public LScore(String[] strArray, double targetmb, String name, double subscore, double penalty) throws IOException {
		alignmentname = name;
		mm = subscore;
		pt = penalty;
		strings = new String[strArray.length];
		factorials = new int[strArray.length];
		factorials[0] = 1;
		for (int i=0; i<strArray.length; i++) strings[i]=","+strArray[i];
		modStrings(targetmb);
		maxmatrix = new double[totalpartnum];
		maxindexmatrix = new int[totalpartnum];
		scoreAll();
		
		BufferedWriter fw = new BufferedWriter(new FileWriter("C:/users/coka/desktop/programmierung/java/projekte/alignments/"+alignmentname+"/parameters.txt"));
		String plty = String.valueOf(pt);
		String mtch = String.valueOf(mt);
		String sub = String.valueOf(mm);
		fw.write("Matchscore: "+mtch+"; Subscore: "+sub+"; Penalty (scalar): "+plty);
		fw.close();
		fw = new BufferedWriter(new FileWriter("C:/users/coka/desktop/programmierung/java/projekte/alignments/"+alignmentname+"/strings.txt"));
		String allst = "";
		for (int i=0; i<strings.length; i++) allst+="str"+String.valueOf(i+1)+":\n"+strings[i]+"\n";
		fw.write(allst);
		fw.close();
		
		
		DataMaster infobj = new DataMaster(alignmentname, strings, partdims, partsize, maxmatrix, maxindexmatrix);
		infobj.pt = pt; infobj.mm = mm; infobj.mt = mt; infobj.segfactorials = factorials;
		String path2 = "C:/users/coka/desktop/programmierung/java/projekte/alignments/"+alignmentname+"/datamaster.object";
		FileOutputStream fos2 = new FileOutputStream(path2);
		ObjectOutputStream oos2 = new ObjectOutputStream(fos2);
		oos2.writeObject(infobj);
		oos2.close();
		
	}
	
	public void scoreAll() {
		// score all segments
		
		Scoring[] partlist;
		int maxdim = 0;
		for (int i=1 ; i<partdims.length; i++) if (partdims[maxdim]<partdims[i]) maxdim=i;
		int layerlength = factorials[strings.length-1]*partdims[strings.length-1]/partdims[maxdim];
		partlist = new Scoring[layerlength];
		int layerstep = 1;
		for (int i=0; i<maxdim; i++) layerstep*=partdims[i];
		int[] subfactorials = factorials.clone();
		for (int i=1; i<factorials.length; i++) 
			if (i>=maxdim) subfactorials[i] = factorials[i]/partdims[maxdim];
		
		
		// schleife, um alle segmente in richtung maxdim abzuarbeiten.. 
		for (int i=0; i<partdims[maxdim]; i++) {
			for (int k=0; k<layerlength; k++) {
				int n = i*layerstep;
				for (int j=0; j<partdims.length; j++) if (j!=maxdim) {
					n += ((k/subfactorials[j])%partdims[j])*factorials[j];
				}
				String[] partStrings = new String[strings.length];
				boolean[] border = new boolean[strings.length];
				for (int j=0; j<strings.length; j++) 
					if ((n/factorials[j])%partdims[j]==partdims[j]-1) {
						partStrings[j] = strings[j].substring(((n/factorials[j])%partdims[j])*partsize[j], ((n/factorials[j])%partdims[j]+1)*(partsize[j]));
						border[j] = true;
					} else {
						partStrings[j] = strings[j].substring(((n/factorials[j])%partdims[j])*partsize[j], ((n/factorials[j])%partdims[j]+1)*(partsize[j])+1);
					}
				Scoring S = new Scoring(partStrings, n, border, alignmentname);
				S.local = true;
				S.pt = pt; S.mt = mt; S.mm = mm;
				
				//System.out.println();
				//System.out.println("i: "+String.valueOf(i)+", k: "+String.valueOf(k)+", n: "+String.valueOf(n));
				
				for (int j=0; j<strings.length; j++)
					if ((n/factorials[j])%partdims[j]>0) {
						int m = n-factorials[j] ,q=0;
						for (int p=0; p<partdims.length; p++) 
							if (p!=maxdim) q+=((m/factorials[p])%partdims[p])*subfactorials[p];
						//System.out.println();
						//System.out.println(n);
						//System.out.println(m);
						//System.out.println(q);
						//if (i==1) System.exit(0);
						S.setFirst(j, partlist[q].borderlist[j]);
					}
				S.score(S.leavef);
				try {
					S.flush();
				} catch (IOException e) {
					e.printStackTrace();
					System.out.println("Problem with flushing");
					System.exit(0);
				}
				maxmatrix[n] = S.globalmax;
				maxindexmatrix[n] = S.globalmaxindex;
				partlist[k] = S;
			}
			System.out.println("-------------------------------------------\nFinished layer "+String.valueOf(i));
		}
		System.out.println("-------------------------------------------");
	}
	
	public void modStrings(double targetmb) {
		int[] y = new int[strings.length];
		double allprod = 1; for (int i=0; i<strings.length; y[i]=strings[i].length(), i++) allprod*=strings[i].length();
		double allmb = allprod*(4+y.length*4)/1000000*1.04; 	//stimmt jetzt grade nur für 2d matrizen, man muss die 4 in der klammer anpassen. (ZEILE 191!)
		double targetparts = allmb/targetmb+1;
		ArrayList<Integer>[] divisors = new ArrayList[y.length];
		int j = 1;
		int[] appendages = new int[strings.length];
		
		for (int n=0; n<y.length; n++) { 
			divisors[n] = new ArrayList<Integer>();
			for (int i=1; i<y[n]/Math.log(y[n]); i++) {
				if (y[n]%i==0) divisors[n].add(i);
			}
			int µ = 0;
			while (divisors[n].size()<10 && strings[n].length()>100) {
				y[n]++;
				String appendix = String.valueOf(Character.toLowerCase(n+µ));
				µ += y.length;
				strings[n] += appendix;
				appendages[n] += 1;
				divisors[n] = new ArrayList<Integer>();
				for (int i=2; i<y[n]/Math.log(y[n]); i++) {
					if (y[n]%i==0) divisors[n].add(i);
				}
			}
			
			if (n>0) factorials[n] = j;
			j *= divisors[n].size(); 
		}
		
		double[] productlogs = new double[j];
		double[] mbdeltas=new double[j], ptsizedeltas=new double[j];
		
		for (int i=0; i<j; i++) {
			long n = 1;
			for (int k=0; k<divisors.length; k++)
				n *= (int)divisors[k].get((int)((double)i/(double)factorials[k])%divisors[k].size());
			double[] partlens = new double[y.length];
			for (int k=0; k<partlens.length; k++) partlens[k] = y[k]/divisors[k].get((i/factorials[k])%divisors[k].size());
			double deltaptsize = 1;
			double minptsize=partlens[0], maxptsize=partlens[0];
			for (int k=1; k<partlens.length; k++) {
				if (minptsize>partlens[k]) minptsize = partlens[k];
				if (maxptsize<partlens[k]) maxptsize = partlens[k];
			} deltaptsize = maxptsize/minptsize;
			
			double deltamb;
			if (n>targetparts) deltamb = n/targetparts;
			else deltamb = targetparts/n;
			
			mbdeltas[i] = deltamb;
			ptsizedeltas[i] = deltaptsize;
		}
		
		for (int i=0; i<mbdeltas.length; i++)
			productlogs[i] = -((Math.pow(mbdeltas[i]-0.5, 2)+0.75)*(ptsizedeltas[i]*ptsizedeltas[i]));
		
		int maxindex = 0;
		for (int i=1; i<productlogs.length; i++) 
			if (productlogs[i]>productlogs[maxindex]) maxindex = i;
		
		partdims = new int[y.length];
		for (int i=0; i<partdims.length; i++) partdims[i] = divisors[i].get((maxindex/factorials[i])%divisors[i].size());
		System.out.print("\nSplit Matrix Dimensions:\n(");
		for (int i=0; i<partdims.length; i++) 
			if (i<partdims.length-1) System.out.print(String.valueOf(partdims[i])+", ");
			else System.out.print(String.valueOf(partdims[i])+")\ninto parts of size\n(");
		partsize = new int[y.length];
		for (int i=0; i<partdims.length; i++) partsize[i] = (int)y[i]/partdims[i];
		for (int i=0; i<partsize.length; i++) 
			if (i<partsize.length-1) System.out.print(String.valueOf(partsize[i])+", ");
			else System.out.print(String.valueOf(partsize[i])+")\n");
		for (int i=0; i<partsize.length; totalpartsize*=partsize[i],totalpartnum*=partdims[i], i++);
		System.out.println("That amounts to "+String.valueOf((double)(totalpartsize)/**...*/*(4+y.length*4)/1000000*1.04)+"MB per part");
		
		System.out.print("Appendages per String: ");
		for (int i=0; i<appendages.length; i++) System.out.print(String.valueOf(appendages[i])+", ");
		System.out.print("\n\n");
		
		for (int i=1; i<partdims.length ; i++) factorials[i]=factorials[i-1]*partdims[i-1];
	}
}
