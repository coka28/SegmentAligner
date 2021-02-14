package Alignment;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;

public class DataMaster implements Serializable {
	
	private static final long serialVersionUID = 1L;
	private ArrayList<Integer> live = new ArrayList<Integer>(); 
	private ArrayList<double[]> scorescache = new ArrayList<double[]>();
	private ArrayList<boolean[][]> pointerscache = new ArrayList<boolean[][]>();
	public ArrayList<boolean[]> excludelistcache = new ArrayList<boolean[]>();  // liste von globalen indices, die ausgelassen werden sollen beim rescoring
	
	
	public double maxmb = 2000;
	public int[] partsize;
	public int[] partdims;
	int[] globalsize;
	public int[] segfactorials;
	public int[] partfactorials;
	public long[] globalfactorials;
	public long[] accesscounter;
	public long partsizeprod;
	public String[] strings;
	double[] maxmatrix;
	int[] maxindexmatrix;
	public double pt,mm,mt;
	
	public String dataname;
	
	public DataMaster(String name, String[] strArray, int[] segdims, int[] segsize, double[] maxmatr, int[] maxindexmatr) {
		maxmatrix = maxmatr;
		maxindexmatrix = maxindexmatr;
		dataname = name;
		partdims = segdims;
		partsize = segsize;
		partfactorials = new int[segdims.length];
		globalfactorials = new long[segdims.length];
		globalsize = new int[segdims.length];
		partfactorials[0] = 1;
		globalfactorials[0] = 1;
		for (int i=1; i<partsize.length; i++) partfactorials[i]=partfactorials[i-1]*partsize[i-1];
		strings = strArray;
		for (int i=1; i<strings.length; i++) globalfactorials[i]=globalfactorials[i-1]*strings[i-1].length();
		partsizeprod = 1;
		int o = 1;
		for (int i=0; i<segdims.length; i++) o*=segdims[i];
		accesscounter = new long[o];
		for (int i=0; i<partsize.length; globalsize[i] = strings[i].length(), i++) partsizeprod*=partsize[i];
	}
	
	public long getMaxIndex(double minscore) throws ClassNotFoundException, IOException {
		long res;
		int maxpos = 0;
		boolean legitmin = false;
		for (int i=0; i<maxmatrix.length; i++) {
			if (maxmatrix[maxpos]<=maxmatrix[i] && maxmatrix[i]>=minscore) {maxpos = i; legitmin = true;}
		}
		
		if (!legitmin) return 0;
		res = maxpos * partsizeprod;
		res += maxindexmatrix[maxpos];
		return res;
	}
	
	public Object[] get(long globalindex, boolean reeval_max) throws ClassNotFoundException, IOException {
		int segmentindex = (int) (globalindex/partsizeprod);
		int localindex = (int) (globalindex%partsizeprod);
		acknowledge(segmentindex, reeval_max);
		
		int tmpindex = live.indexOf(segmentindex);
		Object[] res = {(scorescache.get(tmpindex))[localindex], (pointerscache.get(tmpindex))[localindex], (excludelistcache.get(tmpindex))[localindex]};
		return res;
	}
	
	public void set(long globalindex, double score, boolean[] pointer) throws ClassNotFoundException, IOException {
		int segmentindex = (int) (globalindex/partsizeprod);
		int localindex = (int) (globalindex%partsizeprod);
		acknowledge(segmentindex, true);
		
		int tmpindex = live.indexOf(segmentindex);
		scorescache.get(tmpindex)[localindex] = score;
		pointerscache.get(tmpindex)[localindex] = pointer;
	}

	public double match(long[] position) {
		double res = 0;
		for (int i=0; i<strings.length; i++)
			for (int k=i+1; k<strings.length; k++) {
				if (strings[i].charAt((int)position[i])==strings[k].charAt((int)position[k]))
					res += mt; 
				else res += mm; }
		return res;
	}
	
	boolean[][] legitPointers(long[] position) {
		int[] maxPos = new int[position.length];
		boolean[][] res;
		int n=0, m=1;
		for (int i=0; i<maxPos.length; i++) if (position[i]>0) {
			maxPos[n] = i;
			n++;
			m *= 2;}
		m--;
		res = new boolean[m][maxPos.length];
		int j = 1;
		for (int i=0; i<m; i++,j++) for (int k=0; k<n; k++)
			res[i][maxPos[k]] = ((int)(j/Math.pow(2, k))%2==1);
		//System.out.print("Pointers : ");
		//for (int i=0; i<res.length; i++) print(res[i]);
		return res;
	}
	
	public void boxrescore(int start, int end) throws ClassNotFoundException, IOException {
		int[] chunkdims=new int[segfactorials.length], pos=new int[segfactorials.length], 
				endpos=new int[segfactorials.length], chunkfactorials=new int[segfactorials.length+1];
		int chunksize = 1;
		chunkfactorials[0] = 1;
		for (int i=0; i<segfactorials.length; i++) {
			pos[i] = (start/segfactorials[i])%partdims[i];
			endpos[i] = (end/segfactorials[i])%partdims[i];
			chunkdims[i] = ((end/segfactorials[i])%partdims[i])+1 - ((start/segfactorials[i])%partdims[i]);
			chunksize *= chunkdims[i];
			chunkfactorials[i+1] = chunkfactorials[i] * chunkdims[i];
		}
				
		int[] rescorelist = new int[chunksize]; 
		
		int currentindex = start;
		
		rescorelist[0] = start;
		for (int i=1; i<chunksize; i++) {
			currentindex = start;
			for (int k=0; k<pos.length; k++) currentindex+=((i/chunkfactorials[k])%chunkdims[k])*segfactorials[k];
			
			rescorelist[i] = currentindex;
			
		}
		if (currentindex!=end) {
			System.out.println("ERROR! Rescoring is faulty!");
			System.exit(0);
		}
		
		
		
		System.out.print("\n--rescoring "+String.valueOf(rescorelist.length)+" parts ");
		long t0 = System.currentTimeMillis();
		int counter = 0, progression=1;
		for (int current : rescorelist) {
			acknowledge(current);
			for (int i=0; i<partsizeprod; i++) 
				if (excludelistcache.get(live.indexOf(current))[i]) {
					
				} else {
					long currentglobalindex = globalIndex(current,i);
					long[] currentglobalposition = globalPosition(currentglobalindex);
					boolean[][] p = legitPointers(currentglobalposition);
					
					boolean[] maxPointer = new boolean[currentglobalposition.length];
					double mtc=match(currentglobalposition), max=0;
					double[] scores = new double[p.length];
					for (int k=0; k<p.length; k++) {
						int zeros = 0;
						for (int n=0; n<currentglobalposition.length; n++)
							if (!p[k][n]) zeros++;
						double plty = zeros*pt;
						long newindex = currentglobalindex;
						for (int j=0; j<p[k].length; j++) if (p[k][j]) newindex-=partfactorials[j];
						Object[] tmp = get(newindex);
						double tmpscore = (double)tmp[0];
						scores[k] = tmpscore+mtc-plty;
						if (scores[k]>=max) {
							max = scores[k];
							maxPointer = p[k];
						}
					}
					if (max<=0)	set(currentglobalindex,0,null);
					else set(currentglobalindex,max,maxPointer);
				}
			
			double localmax = 0;
			int maxindex = 0;
			double[] tmp = scorescache.get(live.indexOf(current));
			for (int i=tmp.length-1; i>=0; i--) if (tmp[i]>localmax) {
				localmax = tmp[i];
				maxindex = i;
			}
			
			maxmatrix[current] = localmax;
			maxindexmatrix[current] = maxindex;
			
			counter++;
			double step = ((double)rescorelist.length+0.5)/Math.min(20, rescorelist.length);
			
			if (counter==(int)((double)progression*step)) {
				progression++;
				System.out.print(".");
			}
		}
		System.out.print(" "+String.valueOf((double)(System.currentTimeMillis()-t0)/1000)+"s--\n");
		dumpAll();

	}
		
	private boolean rescore(int segindex, boolean startpoint, boolean force_rescore) throws ClassNotFoundException, IOException { 
		// returns true if all the newly scored entries on the lower borders are identical to the entries in scoringcache,
		// else returns false and rescores specified segment
		acknowledge(segindex, false);
		for (int i=0; i<strings.length; i++) if (segindex>=segfactorials[i]) acknowledge(segindex-segfactorials[i], false);
			
		int cacheindex = live.indexOf(segindex);
		
		double[] scoring = scorescache.get(cacheindex);
		boolean[][] pointers = pointerscache.get(cacheindex);
		boolean[] excludelist = excludelistcache.get(cacheindex);
		boolean[] leavefirst = new boolean[strings.length];
		boolean allthesametop = true;
		boolean allthesamebot = true;
		double localmax = 0;
		int localmaxindex = 0;
		
		
		if (segindex>0) {
			for (int i=0; i<strings.length; i++) 
				if ((segindex/segfactorials[i])%partdims[i]>0) leavefirst[i] = true;
			for (int i=0; i<scoring.length; i++) 
				for (int dim=0; dim<strings.length; dim++)
					if (!startpoint) {
						if ((int)(i/partfactorials[dim])%partsize[dim]==0) 
							if (!excludelist[i]) {
								//double tmp = score(startindex+i, false);
								long globalindex = globalIndex(segindex,i);
								long[] globalpos = globalPosition(globalindex);
								double max = 0;
								double matchscore = match(globalpos);
								boolean[] pointer = new boolean[strings.length];
								boolean[][] p = legitPointers(globalpos);
								for (int k=0; k<p.length; k++) {
									long[] nextpos = globalpos.clone();
									int zeros = 0;
									for (int n=0; n<strings.length; n++)
										if (p[k][n]) nextpos[n]--;
										else zeros++;
									long nextindex = globalIndex(nextpos);
									double score1 = (double)get(nextindex, false)[0];
									double score2 = matchscore-pt*zeros;
									if (score1+score2>max) {
										max = score1+score2;
										pointer = p[k];
									}
									
								}
								if (Math.abs(max-scoring[i])>0.001) 
									allthesametop=false;
								
								if (max>0) {
									scoring[i] = max;
									pointers[i] = pointer;
								} else {
									scoring[i] = 0;
									pointers[i] = null;
								}
								if (localmax<max) {
									localmax = max;
									localmaxindex = i;
								}
								
								
							}
					}
		}
		
		
		double min=0;
		for (int i=0; (!allthesametop || startpoint || force_rescore) && i<scoring.length; i++) {
			boolean skipthis = false;
			for (int k=0; k<leavefirst.length; k++) if (leavefirst[k]
				&& (i/partfactorials[k])%partsize[k]==0) skipthis=true;
			if (!excludelist[i] && !skipthis) {
				long[] currentglobalpos = globalPosition(globalIndex(segindex,i));
				boolean[][] p = legitPointers(currentglobalpos);
				boolean[] maxPointer = new boolean[strings.length];
				double mtc=match(currentglobalpos), max=min;
				int maxindex = 0;
				double[] scores = new double[p.length];
				for (int k=0; k<p.length; k++) {
					int zeros = 0;
					for (int n=0; n<strings.length; n++)
						if (!p[k][n]) zeros++;
					double plty = zeros*pt;
					int nextindex = i;
					for (int n=0; n<strings.length; n++) if (p[k][n]) nextindex-=partfactorials[n];
					scores[k] = scoring[nextindex]+mtc-plty;
					if (scores[k]>=scores[maxindex]) {
						max = scores[k];
						maxPointer = p[k];
						maxindex = k;
					} else if (scores[k]<min) min = scores[k];
				}
				
				for (int dim=0; startpoint && dim<strings.length; dim++) 
					if ((int)(i/partfactorials[dim])%partsize[dim]==partsize[dim]-1)
						if (Math.abs(Math.max(max, 0)-scoring[i])>0.001)
							allthesamebot=false;
				
				if (max<=0) {
					scoring[i] = 0;
					pointers[i] = null;
				} else {
					scoring[i] = max;
					pointers[i] = maxPointer;
				}
			}
			
			if (scoring[i]>localmax) {
				localmax = scoring[i];
				localmaxindex = i;
			}
		}
		
		maxmatrix[segindex] = localmax;
		maxindexmatrix[segindex] = localmaxindex;
		
		cacheindex = live.indexOf(segindex);
		scorescache.set(cacheindex, scoring);
		pointerscache.set(cacheindex, pointers);
		return allthesametop&&allthesamebot;
	}
	
	public void selective_rescore(ArrayList<Integer> startlist) throws ClassNotFoundException, IOException {
		ArrayList<Integer> rescorelist = new ArrayList<Integer>();
		
		rescorelist.addAll(startlist);
		for (int i=0; i<rescorelist.size(); i++)
			for (int k=0; k<rescorelist.size()-1; k++) 
				if (rescorelist.get(k)>rescorelist.get(k+1)) {
					int tmp = rescorelist.get(k);
					rescorelist.set(k, rescorelist.get(k+1));
					rescorelist.set(k+1, tmp);
				}

		System.out.print("-------------------------------------------\n--rescoring: "+String.valueOf(rescorelist.get(0))+", ");
		boolean morethanone = !rescore(rescorelist.get(0),true,true);
		for (int i=0; i<strings.length; i++)
			if (morethanone && (rescorelist.get(0)/segfactorials[i])%partdims[i]<partdims[i]-1 
					&& !rescorelist.contains(rescorelist.get(0)+segfactorials[i])) 
				rescorelist.add(rescorelist.get(0)+segfactorials[i]);
		
		rescorelist.remove(0);
		boolean force;
		
		while (rescorelist.size()>0) {
			for (int i=0; i<rescorelist.size(); i++)
				for (int k=0; k<rescorelist.size()-1; k++) 
					if (rescorelist.get(k)>rescorelist.get(k+1)) {
						int tmp = rescorelist.get(k);
						rescorelist.set(k, rescorelist.get(k+1));
						rescorelist.set(k+1, tmp);
					}
			
			System.out.print(String.valueOf(rescorelist.get(0))+", ");

			int current = rescorelist.get(0);
			force = startlist.contains(current);
			
			if (!rescore(current, false, force)) {
				for (int i=0; i<strings.length; i++)
					if ((current/segfactorials[i])%partdims[i]<partdims[i]-1 && !rescorelist.contains(current+segfactorials[i])) 
						rescorelist.add(current+segfactorials[i]);
			}
			
			rescorelist.remove(rescorelist.indexOf(current));
			
			
		}
		System.out.println("finished--\n-------------------------------------------");
		dumpAll(false);
	}

	public long globalIndex(int segindex, int localindex) {
		return (long)(segindex)*partsizeprod+localindex;
	}
	
	public long globalIndex(long[] globalposition) {
		int[] localpos = new int[globalposition.length];
		for (int i=0; i<globalposition.length; i++) localpos[i]=(int)(globalposition[i]%partsize[i]);
		int localindex = 0;
		for (int i=0; i<globalposition.length; i++) localindex+=localpos[i]*partfactorials[i];
		return globalIndex(segmentIndex(globalposition), localindex);
	}
	
	public long[] globalPosition(long globalindex) {
		long[] res = new long[strings.length];
		int seg = segmentIndex(globalindex);
		for (int i=0; i<res.length; i++) res[i]=((seg/segfactorials[i])%partdims[i])*partsize[i];
		int locindex = (int)(globalindex%partsizeprod);
		for (int i=0; i<res.length; i++) res[i]+=(locindex/partfactorials[i])%partsize[i];
		return res;
	}
	
	public int segmentIndex(long globalindex) {
		return (int)(globalindex/partsizeprod);
	}
	
	public int segmentIndex(long[] globalposition) {
		int res = 0;
		int[] segpos = new int[globalposition.length];
		for (int i=0; i<globalposition.length; i++) segpos[i]=(int)(globalposition[i]/partsize[i]);
		for (int i=0; i<globalposition.length; i++) res+=segpos[i]*segfactorials[i];
		return res;
	}
	
	public void erase(long globalindex) throws ClassNotFoundException, IOException {
		int segmentindex = (int) (globalindex/partsizeprod);
		int localindex = (int) (globalindex%partsizeprod);
		acknowledge(segmentindex, true);
		
		int tmpindex = live.indexOf(segmentindex);
		scorescache.get(tmpindex)[localindex] = 0;
		pointerscache.get(tmpindex)[localindex] = null;
		excludelistcache.get(tmpindex)[localindex] = true;
	}
	
	private void acknowledge(int segmentindex, boolean reeval_max) throws ClassNotFoundException, IOException {
		if (!live.contains(segmentindex)) {
			loadSegment(segmentindex);
		} else if (!(live.indexOf(segmentindex)<live.size()-2)) {
			int tmp = live.indexOf(segmentindex);
			int segindextmp = live.get(tmp);
			double[] scoringtmp = scorescache.get(tmp);
			boolean[][] pointerstmp = pointerscache.get(tmp);
			boolean[] excludetmp = excludelistcache.get(tmp);
			live.remove(tmp);
			scorescache.remove(tmp);
			pointerscache.remove(tmp);
			excludelistcache.remove(tmp);
			live.add(segindextmp);
			scorescache.add(scoringtmp);
			pointerscache.add(pointerstmp);
			excludelistcache.add(excludetmp);
		}
		
		accesscounter[segmentindex] = -1;
		while (live.size()>=maxmb/(1.25e-5*partsizeprod)) {
			if (live.size()<=3) {
				System.out.println("Can simultaneously hold only 3 or less segments with maxmb capacity: exiting for data safety");
				System.exit(0);
			}
			dumpSegment(live.get(0), reeval_max);
		}
		
	}
	
	public void dumpAll(boolean reeval_max) throws IOException {
		while (live.size()>0) {
			int lastindex = live.size()-1;
			int index = live.get(lastindex);
			if (reeval_max) {
				double max = 0;
			    int maxindex = 0;
			    double[] scores = scorescache.get(lastindex);
			    for (int i=0; i<scores.length; i++) 
			    	if (scores[i]>max) {
			    		max = scores[i];
			    		maxindex = i;
			    	}
			    maxmatrix[index] = max;
			    maxindexmatrix[index] = maxindex;
			}
			String fileName = "C:/users/coka/desktop/programmierung/java/projekte/alignments/"+dataname + "/" + String.valueOf(live.get(lastindex));
		    FileOutputStream fos = new FileOutputStream(fileName);
		    ObjectOutputStream oos = new ObjectOutputStream(fos);
		    Object[] objtmp = {scorescache.get(lastindex), pointerscache.get(lastindex), excludelistcache.get(lastindex)};
		    oos.writeObject(objtmp);
		    live.remove(lastindex);
		    scorescache.remove(lastindex);
		    pointerscache.remove(lastindex);
		    excludelistcache.remove(lastindex);
		    oos.close();
		}
	}
	
	private void dumpSegment(int index, boolean reeval_max) throws IOException, ClassNotFoundException {
	    String fileName = "C:/users/coka/desktop/programmierung/java/projekte/alignments/"+dataname + "/" + String.valueOf(index);
	    FileOutputStream fos = new FileOutputStream(fileName);
	    ObjectOutputStream oos = new ObjectOutputStream(fos);
	    int tmp = live.indexOf(index);
	    if (reeval_max) {
		    double max = 0;
		    int maxindex = 0;
		    double[] scores = scorescache.get(tmp);
		    for (int i=0; i<scores.length; i++) 
		    	if (scores[i]>max) {
		    		max = scores[i];
		    		maxindex = i;
		    	}
		    maxmatrix[index] = max;
		    maxindexmatrix[index] = maxindex;
	    }
	    Object[] objtmp = {scorescache.get(tmp), pointerscache.get(tmp), excludelistcache.get(tmp)};
	    oos.writeObject(objtmp);
	    live.remove(tmp);
	    scorescache.remove(tmp);
	    pointerscache.remove(tmp);
	    excludelistcache.remove(tmp);
	    oos.close();
	}
	
	private void loadSegment(int index) throws IOException, ClassNotFoundException {
		String fileName = "C:/users/coka/desktop/programmierung/java/projekte/alignments/"+dataname + "/" + String.valueOf(index);
		FileInputStream fin = new FileInputStream(fileName);
		ObjectInputStream ois = new ObjectInputStream(fin);
		live.add(index);
		Object[] tmp = (Object[]) ois.readObject();
		scorescache.add((double[])(tmp[0]));
		pointerscache.add((boolean[][])(tmp[1]));
		excludelistcache.add((boolean[])(tmp[2]));
		ois.close();
	}

	public void reMax(ArrayList<Integer> list) throws ClassNotFoundException, IOException {
		for (int i: list) {
			acknowledge(i,false);
			double[] scores = scorescache.get(live.indexOf(i));
			double max = 0; int maxindex = 0;
			for (int k=0; k<scores.length; k++) 
				if (max<scores[k]) {
					max = scores[k];
					maxindex = k;
				}
			
			maxmatrix[i] = max;
			maxindexmatrix[i] = maxindex;
		}
	}
	
}
