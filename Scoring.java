package Alignment;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;

public class Scoring {

	public double pt=2, mt=1, mm=-1;
	public double[] scoring;
	public boolean[][] pointers;
	public String[] strings;
	public int[] strLengths, factorials;
	public boolean local = false;
	public int globalmaxindex=1, maxIndex=1;
	public double globalmax=0;
	public boolean[] borders;
	public boolean[] leavef;
	public int segmentindex = 0;
	public Object[][] borderlist;
	String dataname;
	public int[] excludelist;
	
	public Scoring(String[] args) {
		strings = new String[args.length];
		strLengths = new int[args.length];
		factorials = new int[args.length];
		factorials[0] = 1;
		for (int i=0; i<strings.length; i++) {
			strings[i] = "," + args[i];
			strLengths[i] = strings[i].length();
			maxIndex *= strLengths[i];
			if (i>0) factorials[i] = factorials[i-1]*strLengths[i-1];
		}
		
		scoring = new double[maxIndex];
		pointers = new boolean[maxIndex][];
	}
	
	public Scoring(String[] args, int segindex, boolean[] border, String name) {
		borders = border;
		dataname = name;
		segmentindex = segindex;
		leavef = new boolean[args.length];
		borderlist = new Object[args.length][];
		strings = new String[args.length];
		strLengths = new int[args.length];
		factorials = new int[args.length];
		factorials[0] = 1;
		for (int i=0; i<strings.length; i++) {
			strings[i] = args[i];
			strLengths[i] = strings[i].length();
			maxIndex *= strLengths[i];
			if (i>0) factorials[i] = factorials[i-1]*strLengths[i-1];
		}
		
		scoring = new double[maxIndex];
		pointers = new boolean[maxIndex][];
	}
	
	public void setFirst(int dim, Object[] borders) {
		leavef[dim] = true;
		
		double[] scoringborder = (double[]) borders[0];
		boolean[][] pointerborder = (boolean[][]) borders[1];
		for (int i=0,k=0; i<scoring.length; i++) {
			if ((int)(i/factorials[dim])%strLengths[dim]==0) {
				scoring[i] = scoringborder[k];
				pointers[i] = pointerborder[k];
				k ++;
			}
		}
	}
	
	public void score(boolean[] leavefirst) {		
		double min=0;
		for (int i=0; i<maxIndex; i++) {
			boolean skipthis = false;
			for (int k=0; k<leavefirst.length; k++) if (leavefirst[k]
				&& (i/factorials[k])%strLengths[k]==0) skipthis=true;
			if (!skipthis) {
				int[] thisPos = position(i);
				boolean[][] p = legitPointers(thisPos);
				boolean[] maxPointer = new boolean[thisPos.length];
				double mtc=match(thisPos), max=min;
				int maxindex = 0;
				double[] scores = new double[p.length];
				for (int k=0; k<p.length; k++) {
					int zeros = 0;
					for (int n=0; n<thisPos.length; n++)
						if (!p[k][n]) zeros++;
					double plty = zeros*pt;
					scores[k] = scoring[index(step(thisPos,p[k]))]+mtc-plty;
					if (scores[k]>=scores[maxindex]) {
						max = scores[k];
						maxPointer = p[k];
						maxindex = k;
					} else if (scores[k]<min) min = scores[k];
				}
				if (local && max<=0) {
					scoring[i] = 0;
					pointers[i] = null;
				} else {
					scoring[i] = max;
					pointers[i] = maxPointer;
				}
			}
			
			if (scoring[i]>globalmax) {
				globalmax = scoring[i];
				globalmaxindex = i;
			}
		}
	}

	Object[] getLast(int dim) { 	//return scoring and pointer matrix border
		Object[] res = new Object[2];
		int bdlength = 1;
		for (int i=0; i<strLengths.length; i++) if (i!=dim) bdlength*=strLengths[i];
		double[] scoreborder = new double[bdlength];
		boolean[][] pointerborder = new boolean[bdlength][];
		 
		for (int i=0,k=0; i<scoring.length; i++) {
			if ((int)(i/factorials[dim])%strLengths[dim]==strLengths[dim]-1) {
				scoreborder[k] = scoring[i];
				pointerborder[k] = pointers[i];
				k ++;
			}
		}
		res[0] = scoreborder;
		res[1] = pointerborder;
		return res;
	}
	
	public void flush() throws IOException { // to cut last indices and save matrices in file and void local data, except for borders
		int[] dims;
		int a = 0;
		for (int i=0; i<borders.length; i++) if (!borders[i]) {
			borderlist[i] = getLast(i);
			a++;
		}
		dims = new int[a];
		for (int i=0,k=0; i<borders.length; i++) if (!borders[i]) {
			dims[k] = i;
			k++;
		}
		
		int j = 0;
		for (int i=0; i<scoring.length; i++) {
			boolean delete = false;
			for (int k=0; k<dims.length; k++) 
				if ((i/factorials[dims[k]])%strLengths[dims[k]]==strLengths[dims[k]]-1)
					delete = true;
			if (!delete) j++;
		}
		globalmax = 0;
		globalmaxindex = 0;
		double[] newScoring = new double[j];
		boolean[][] newPointers = new boolean[j][];
		
		j = 0;
		for (int i=0; i<scoring.length; i++) {
			boolean insert = true;
			for (int k=0; k<dims.length; k++) if ((i/factorials[dims[k]])%strLengths[dims[k]]==strLengths[dims[k]]-1) insert = false;
			if (insert) {
				newScoring[j] = scoring[i];
				newPointers[j] = pointers[i];
				j++;
			}
		}
		
		for (int i=newScoring.length-1; i>=0; i--) if (globalmax<newScoring[i]) {
			globalmax = newScoring[i];
			globalmaxindex = i;
		}
		
		for (int i=0; i<dims.length; i++) {
			strLengths[dims[i]]--;
			strings[dims[i]] = strings[dims[i]].substring(0, strings[dims[i]].length()-1);
			
		}
		
		for (int i=1; i<strLengths.length; i++)
			factorials[i] = factorials[i-1] * strLengths[i-1];
		
		boolean[] excludematrix = new boolean[newScoring.length];
		
		Object[] packaged = {newScoring, newPointers, excludematrix};
		String fileName = "C:/users/coka/desktop/programmierung/java/projekte/alignments/"+dataname+"/"+String.valueOf(segmentindex);
		FileOutputStream fos = new FileOutputStream(fileName);
		ObjectOutputStream oos = new ObjectOutputStream(fos);
		oos.writeObject(packaged);
		oos.close();
		
		scoring = null;
		pointers = null; 
		newScoring = null;
		newPointers = null;
		System.gc();
	}

	double match(int[] position) {
		double res = 0;
		for (int i=0; i<strings.length; i++)
			for (int k=i+1; k<strings.length; k++) {
				if (strings[i].charAt(position[i])==strings[k].charAt(position[k]))
					res += mt; 
				else res += mm; }
		return res;
	}
	
	int[] step(int[] a, boolean[] b) {
		int[] res = new int[a.length];
		for (int i=0; i<a.length; i++)
			if (b[i]) res[i] = a[i]-1; 
			else res[i] = a[i];
		return res;
	}
	
	int index(int[] position) {
		int res = 0;
		for (int i=0; i<position.length; i++) {
			int tmp = 1;
			for (int k=0; k<i; tmp*=strLengths[k], k++);
			res += tmp*position[i];
		}
		return res;
	}
	
	int[] position(int index) {
		int[] res = new int[strings.length];
		for (int i=0; i<res.length; i++) {
			int j = 1;
			for (int k=0; k<i; j*=strLengths[k], k++);
			int tmp = index/j;
			res[i] = tmp%strLengths[i];
		}
		return res;
	}
	
	boolean[][] legitPointers(int[] position) {
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

}
