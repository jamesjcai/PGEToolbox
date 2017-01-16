/* this is the control class with the main method. it takes in a file containing all
   samples. it calls different classes to do specific job and output the results to
   different files.*/

import java.util.*;
import java.io.*;
import java.lang.*;

public class EMLD {

	public static void main (String[] args) {
		String curLine;
		int sampleSize;
		int markerNum;
		StringTokenizer St;
		int genotype[][][];
		double D, r2, Dprime;


	try {
		LineNumberReader FR = new LineNumberReader(new FileReader (args[0]));
		curLine = FR.readLine();
		St = new StringTokenizer(curLine);
		sampleSize = Integer.parseInt(St.nextToken());
		markerNum = Integer.parseInt(St.nextToken());
		genotype = new int[sampleSize][markerNum][2];
		int n = 0;
		while ((curLine = FR.readLine()) != null) {
			St = new StringTokenizer(curLine);
			St.nextToken();
			for (int i = 0; i < markerNum; i++) {
				genotype[n][i][0] = Integer.parseInt(St.nextToken());
				genotype[n][i][1] = Integer.parseInt(St.nextToken());
			}
			n++;

		}

		PrintWriter HapFreq = new PrintWriter (new BufferedWriter (new FileWriter("HapFreq.txt")));
		HapFreq.println ("M1\tM2\t11\t12\t21\t22\tsample_size");
		PrintWriter gold = new PrintWriter (new BufferedWriter (new FileWriter("LD.xt")));
		gold.println ("M1\tM2\tD\tD'\tr2");

		for (int i = 0; i< markerNum; i++) {
			for (int j= i+1; j < markerNum; j++) {
				int P1111=0, P1112=0, P1122=0;
				int P1211=0, P1212=0, P1222=0;
				int P2211=0, P2212=0, P2222=0;
				for (int k = 0; k < sampleSize; k++) {
					if (genotype[k][i][0]!=0 && genotype[k][i][1]!=0 &&
						genotype[k][j][0]!=0 && genotype[k][j][1]!=0) {
						if (genotype[k][i][0]==1 && genotype[k][i][1]==1) {
							if (genotype[k][j][0]==1 && genotype[k][j][1]==1) {
								P1111++;
							}
							else if (genotype[k][j][0]==2 && genotype[k][j][1]==2) {
								P1122++;
							}
							else {
								P1112++;
							}
						}
						else if (genotype[k][i][0]==2 && genotype[k][i][1]==2) {
							if (genotype[k][j][0]==1 && genotype[k][j][1]==1) {
								P2211++;
							}
							else if (genotype[k][j][0]==2 && genotype[k][j][1]==2) {
								P2222++;
							}
							else {
								P2212++;
							}
						}
						else {
							if (genotype[k][j][0]==1 && genotype[k][j][1]==1) {
								P1211++;
							}
							else if (genotype[k][j][0]==2 && genotype[k][j][1]==2) {
								P1222++;
							}
							else {
								P1212++;
							}
						}
					}

				}
				EMHaplo EH = new EMHaplo(i, j);
				EH.set(P1111, 1111);
				EH.set(P1112, 1112);
				EH.set(P1122, 1122);
				EH.set(P1211, 1211);
				EH.set(P1212, 1212);
				EH.set(P1222, 1222);
				EH.set(P2211, 2211);
				EH.set(P2212, 2212);
				EH.set(P2222, 2222);
				EH.estimate();

				double P11 = EH.FP11;
				double p1 = EH.p1;
				double p2 = EH.p2;
				double q1 = EH.q1;
				double q2 = EH.q2;

				D = P11 - p1*q1;

				//compute r2;
				r2 = (D*D)/(p1*q1*p2*q2);


				//compute Dprime;
				double Dmax;
				if (D > 0) {
					Dmax = Math.min(p1*q2, p2*q1);
				}else {
					Dmax = Math.min(p1*q1, p2*q2);
				}
				Dprime = D/Dmax;
				gold.println((i+1)+"\t"+(j+1)+"\t"+D+"\t"+Math.abs(Dprime)+"\t"+r2);
				HapFreq.println((i+1)+"\t"+(j+1)+"\t"+EH.FP11+"\t"+EH.FP12+"\t"+EH.FP21+"\t"+EH.FP22+"\t"+EH.n);
				//System.out.println((i+1)+" "+(j+1)+" "+ r2+" "+Math.abs(Dprime));
			}
		}
		gold.close();
		HapFreq.close();
	}catch (Exception e){
		e.printStackTrace();
	}


	}
}


