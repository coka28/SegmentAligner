package App;

import Alignment.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class Exec {

	@SuppressWarnings("unused")
	public static void main(String[] args) throws ClassNotFoundException, IOException {
		
		double target_filesize = 40; 		// in megabytes
		int minimum_score = 5;
		// matchscore is 1
		double subscore = -5;
		double penalty = 5;
		/*
		String[] filepaths = new String[args.length/2];
		String[] stringnames = new String[args.length/2];
		for (int i=0; i<args.length; i+=2) {
			filepaths[i] = args[i];
			stringnames[i] = args[i+1];
		}
		//pairwise
		int n = filepaths.length;
		String[][] pairs = new String[(n*(n-1))/2][];
		String[] titles = new String[(n*(n-1))/2];
		for (int i=0,j=0; i<filepaths.length-1; i++)
			for (int k=i+1; k<filepaths.length; k++,j++) {
				String[] tmp = {filepaths[i], filepaths[k]};
				pairs[j] = tmp;
				String tmp2 = stringnames[i]+"_"+stringnames[k]+"_bis"+String.valueOf(minimum_score);
				titles[j] = tmp2;
			}
		
		for (int i=0; i<pairs.length; i++)
			align(titles[i],pairs[i], target_filesize, minimum_score,subscore,penalty);
		**/
		String sarscov1path = "C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\SARS-CoV.txt"; 
		String sarscov2path = "C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\SARS-CoV-2.txt";
		String merscovpath = "C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\MERS-CoV.txt";
		String HIVpath = "C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\HIV.txt";
		
		String[] paths = {sarscov2path, sarscov1path};
		String[] paths2 = {sarscov2path, merscovpath};
		String[] paths3 = {merscovpath, sarscov1path};
		
		String[] strNames = {"SARS-CoV2","SARS-CoV1","MERS-CoV","HIV1"};
		
		String[] tmp = {"C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\allParts\\sars1_allParts.txt",
				"C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\allParts\\mers_allParts.txt",
				"C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\allParts\\hiv_allParts.txt"};
		//new FindSimilarAlignments("C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\allParts\\sars2_allParts.txt",
				//tmp,strNames);
		
		//align("SARS-CoV-2_SARS-CoV_bis5", paths, target_filesize, minimum_score, subscore, penalty);
		//align("SARS-CoV-2_MERS-CoV_bis5", paths2, target_filesize, minimum_score, subscore, penalty);
		align("MERS-CoV_SARS-CoV_bis5__", paths3, target_filesize, minimum_score, subscore, penalty);
		
		//virusTripleAlignments(target_filesize, minimum_score, subscore, penalty);
		//virusAlignments(target_filesize, minimum_score, subscore, penalty);
		
		//long[] testStringLengths = {342,249,331};
		//globalAlignmentTest(testStringLengths);
		//scoringTimeTest(testStringLengths, false);
		//randomTest(testStringLengths);
		
	}
	
	static void globalAlignmentTest(long[] strLengths) throws IOException {
		for (int k=0; k<strLengths.length; k++) System.out.print(String.valueOf(strLengths[k])+", ");System.out.println();
		String[] str = rand(strLengths);
		long t0 = System.currentTimeMillis();
		new GlobAl(str,1,-1,2);
		System.out.println("Laufzeit: "+String.valueOf(((double)System.currentTimeMillis()-t0)/1000)+"s");
	}
	
	static void scoringTimeTest(long[] strLengths, boolean printmatrix) {
		double t0 = System.currentTimeMillis();
		Scoring a = new Scoring(rand(strLengths));
		a.local = true;
		a.score(null);
		System.out.println((System.currentTimeMillis()-t0)/1000);    
		for (int i=0; printmatrix && i<a.scoring.length; i++) {
			if (i%(strLengths[0]+1)==0) System.out.print(System.lineSeparator()+String.valueOf(a.scoring[i])+"  ");
			else System.out.print(String.valueOf(a.scoring[i])+"  ");
		}
	}
	
	static void randomTest(long[] strLengths) throws ClassNotFoundException, IOException {
		String dataname = "randomtest";
		
		File f = new File("C:/users/coka/desktop/programmierung/java/projekte/alignments/"+dataname);
		if (!f.exists()) f.mkdir();
		else {
			String[]entries = f.list();
			for(String s: entries){
			    File currentFile = new File(f.getPath(),s);
			    currentFile.delete();
			}
		}
		
		long t0 = System.currentTimeMillis();
		for (int k=0; k<strLengths.length; k++) System.out.print(String.valueOf(strLengths[k])+", ");System.out.println();
		
		new LScore(rand(strLengths), 3, dataname, -9, 9);
		System.out.println("Laufzeit für das Scoring: "+String.valueOf((double)(System.currentTimeMillis()-t0)/1000)+"s");
		
		Ltb tb = new Ltb(dataname);
		tb.loadDataMaster();
		tb.getAll(10);
		
		System.out.println("Laufzeit gesamt: "+String.valueOf((double)(System.currentTimeMillis()-t0)/1000)+"s");
	}
	
	static void align(String title, String[] filepaths, double targetsize, double minscore, double subscore, double penalty) throws IOException, ClassNotFoundException {
		File[] files = new File[filepaths.length];
		String[] genomes = new String[filepaths.length];
		for (int i=0; i<filepaths.length; i++) { 
			files[i] = new File(filepaths[i]);
			BufferedReader br = new BufferedReader(new FileReader(files[i]));
			String tmp;
			genomes[i] = "";
			while ((tmp = br.readLine()) != null) genomes[i] += tmp;
			br.close();
		}

		System.out.println(title+"\n");
		File f = new File("C:/users/coka/desktop/programmierung/java/projekte/alignments/"+title);
		if (!f.exists()) f.mkdir();
		else {
			System.out.println("");
			String[] entries = f.list();
			for(String s: entries){
			    File currentFile = new File(f.getPath(),s);
			    currentFile.delete();
			}
		}
		long t0 = System.currentTimeMillis();			
		for (int k=0; k<genomes.length; k++) {
			System.out.print("str"+String.valueOf(k+1)+": "+String.valueOf(genomes[k].length())); 
			System.out.println();
		}
		
		System.out.println();
		for (String st: genomes) System.out.println(st);
		
		new LScore(genomes, targetsize, title, subscore, penalty);
		System.out.println("Laufzeit für das Scoring: "+String.valueOf((double)(System.currentTimeMillis()-t0)/1000)+"s");
		
		Ltb tb = new Ltb(title);
		tb.loadDataMaster();
		tb.getAll(minscore);
		
		System.out.println(title+": Laufzeit: "+String.valueOf((double)(System.currentTimeMillis()-t0)/1000)+"s\n\n");
	}
	
	static void virusTripleAlignments(double filesize, int minScore, double sub, double plty) throws IOException, ClassNotFoundException {
		long starttime = System.currentTimeMillis();
		
		File sars1 = new File("C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\SARS-CoV.txt"); 
		File sars2 = new File("C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\SARS-CoV-2.txt");
		File mers = new File("C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\MERS-CoV.txt");
		
		String tmp;
		BufferedReader br = new BufferedReader(new FileReader(sars1));
		String sars1str = "";
		while ((tmp = br.readLine()) != null) sars1str += tmp;
		br.close();
		
		br = new BufferedReader(new FileReader(sars2));
		String sars2str = "";
		while ((tmp = br.readLine()) != null) sars2str += tmp;
		br.close();
		
		br = new BufferedReader(new FileReader(mers));
		String mersstr = "";
		while ((tmp = br.readLine()) != null) mersstr += tmp;
		br.close();
		
		String[] strings = {sars1str, sars2str, mersstr};
		
		String data_title = "SARS-CoV1_SARS-CoV2_MERS-CoV";
		
		// multiple alignment of SARS-CoV, SARS-CoV2 and MERS-CoV
		System.out.println(data_title+"\n");
		File f = new File("C:/users/coka/desktop/programmierung/java/projekte/alignments/"+data_title);
		if (!f.exists()) f.mkdir();
		else {
			String[] entries = f.list();
			for(String s: entries){
			    File currentFile = new File(f.getPath(),s);
			    currentFile.delete();
			}
		}
			
		long t0 = System.currentTimeMillis();
		for (int k=0; k<strings.length; k++) {
			System.out.print("str"+String.valueOf(k+1)+": "+String.valueOf(strings[k].length())); 
			System.out.println();
		}
		
		System.out.println();
		for (String st: strings) System.out.println(st);
		
		LScore L = new LScore(strings, filesize, data_title, sub, plty);
		System.out.println("Laufzeit für das Scoring: "+String.valueOf((double)(System.currentTimeMillis()-t0)/1000)+"s");
		
		Ltb tb = new Ltb(data_title);
		tb.loadDataMaster();
		tb.getAll(minScore);
		
		System.out.println(data_title+": Laufzeit: "+String.valueOf((double)(System.currentTimeMillis()-t0)/1000)+"s\n\n");
		
		System.out.println("Gesamte Laufzeit für alle Alignments:\n"+String.valueOf((double)(System.currentTimeMillis()-starttime)/1000)+"s\n");
	}
	
	static void virusAlignments(double filesize, int minScore, double sub, double plty) throws IOException, ClassNotFoundException {
		long starttime = System.currentTimeMillis();
		
		File sars1 = new File("C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\SARS-CoV.txt"); 
		File sars2 = new File("C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\SARS-CoV-2.txt");
		File mers = new File("C:\\Users\\Coka\\Desktop\\Programmierung\\Python\\Projekte\\Viruses\\MERS-CoV.txt");
		
		String tmp;
		BufferedReader br = new BufferedReader(new FileReader(sars1));
		String sars1str = "";
		while ((tmp = br.readLine()) != null) sars1str += tmp;
		br.close();
		
		br = new BufferedReader(new FileReader(sars2));
		String sars2str = "";
		while ((tmp = br.readLine()) != null) sars2str += tmp;
		br.close();
		
		br = new BufferedReader(new FileReader(mers));
		String mersstr = "";
		while ((tmp = br.readLine()) != null) mersstr += tmp;
		br.close();
		
		String[] sars_1und2 = {sars1str, sars2str};
		String[] sars1_mers = {sars1str, mersstr};
		String[] sars2_mers = {sars2str, mersstr};
		
		String[] data_titles = {"sars1_sars2", "sars1_mers", "sars2_mers"};
		String[][] data_strings = {sars_1und2, sars1_mers, sars2_mers};
		
		// pairwise alignment of SARS-CoV, SARS-CoV2 and MERS-CoV
		
		for (int i=0; i<data_titles.length; i++) {
			System.out.println(data_titles[i]+"\n");
			File f = new File("C:/users/coka/desktop/programmierung/java/projekte/alignments/"+data_titles[i]);
			if (!f.exists()) f.mkdir();
			else {
				String[] entries = f.list();
				for(String s: entries){
				    File currentFile = new File(f.getPath(),s);
				    currentFile.delete();
				}
			}
			
			long t0 = System.currentTimeMillis();
			for (int k=0; k<data_strings[i].length; k++) {
				System.out.print("str"+String.valueOf(k+1)+": "+String.valueOf(data_strings[i][k].length())); 
				System.out.println();
			}
			
			System.out.println();
			for (String st: data_strings[i]) System.out.println(st);
			
			LScore L = new LScore(data_strings[i], filesize, data_titles[i], sub, plty);
			System.out.println("Laufzeit für das Scoring: "+String.valueOf((double)(System.currentTimeMillis()-t0)/1000)+"s");
			
			Ltb tb = new Ltb(data_titles[i]);
			tb.loadDataMaster();
			tb.getAll(minScore);
			
			System.out.println(data_titles[i]+": Laufzeit: "+String.valueOf((double)(System.currentTimeMillis()-t0)/1000)+"s\n\n");
		}
		
		System.out.println("Gesamte Laufzeit für alle Alignments:\n"+String.valueOf((double)(System.currentTimeMillis()-starttime)/1000)+"s\n");
	}
	
	static String[] rand(long[] lengths) {
		String alphabet = "ATCG";
		String[] res = new String[lengths.length];
		for (int i=0; i<lengths.length; i++) {
			res[i] = "";
			int x = 0;
			for (long k=0; k<lengths[i]; k++) {
				x = (int)(Math.random()*alphabet.length());
				res[i] += alphabet.charAt(x);
			}
		}
		return res;
	}

}
