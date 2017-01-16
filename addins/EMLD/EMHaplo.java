
import java.util.*;
import java.io.*;
import java.lang.*;

public class EMHaplo {

	int P1111=0, P1112=0, P1122=0;
	int P1211=0, P1212=0, P1222=0;
	int P2211=0, P2212=0, P2222=0;
	int locus1, locus2;
	public double p1, p2, q1, q2;
	public double FP11, FP12, FP21, FP22;
	double P11, P12, P21, P22;
	double EP11, EP12, EP21, EP22;
	public int n; //number of individuals

	public EMHaplo (int locus1, int locus2) {
		this.locus1 = locus1;
		this.locus2 = locus2;
	}

	public void estimate() {
		n = P1111+P1112+P1122+P1211+P1212+P1222+P2211+P2212+P2222;
		p1 = (double)(2*(P1111+P1112+P1122)+P1211+P1212+P1222)/(2*n);
		p2 = (double)(2*(P2211+P2212+P2222)+P1211+P1212+P1222)/(2*n);
		q1 = (double)(2*(P1111+P1211+P2211)+P1112+P1212+P2212)/(2*n);
		q2 = (double)(2*(P1122+P1222+P2222)+P1112+P1212+P2212)/(2*n);
		double maxLK = Double.NEGATIVE_INFINITY;
		for (int i = 0; i <100; i++) {
			int maxIt = 1000; //maximum number of iteration;
			double x = Math.random(); //denote the proportion of phenotype 11/22 for double heterozygotes.
			P11 = (double)(2*P1111+P1112+P1211+x*P1212)/(2*n);
			P12 = (double)(P1112+2*P1122+P1222+(1-x)*P1212)/(2*n);
			P21 = (double)(2*P2211+P1211+P2212+(1-x)*P1212)/(2*n);
			P22 = (double)(P1222+P2212+2*P2222+x*P1212)/(2*n);

			double ex = (P11*P22)/(P11*P22 + P12*P21);//E-step
			//M-step;
			EP11 = (double)(2*P1111+P1112+P1211+ex*P1212)/(2*n);
			EP12 = (double)(P1112+2*P1122+P1222+(1-ex)*P1212)/(2*n);
			EP21 = (double)(2*P2211+P1211+P2212+(1-ex)*P1212)/(2*n);
			EP22 = (double)(P1222+P2212+2*P2222+ex*P1212)/(2*n);

			while ((((EP11-P11)>0.000001)||((EP12-P12)>0.000001)
					 ||((EP21-P21)>0.000001)||((EP22-P22)>0.000001)) && maxIt > 0) {
				//System.out.println(EP11+" "+EP12+" "+EP21+" "+EP22);
				ex = (EP11*EP22)/(EP11*EP22+EP12*EP21);
				P11 = EP11;
				P12 = EP12;
				P21 = EP21;
				P22 = EP22;
				EP11 = (double)(2*P1111+P1112+P1211+ex*P1212)/(2*n);
				EP12 = (double)(P1112+2*P1122+P1222+(1-ex)*P1212)/(2*n);
				EP21 = (double)(2*P2211+P1211+P2212+(1-ex)*P1212)/(2*n);
				EP22 = (double)(P1222+P2212+2*P2222+ex*P1212)/(2*n);
				maxIt--;
			}
			double LK = Math.log(Math.pow(P11*P11, P1111))+Math.log(Math.pow(2*P11*P12, P1112))+Math.log(Math.pow(P12*P12, P1122))
			   + Math.log(Math.pow(2*P11*P21, P1211))+Math.log(Math.pow(2*P11*P22+2*P12*P21, P1212))+Math.log(Math.pow(2*P12*P22, P1222))
			   + Math.log(Math.pow(P21*P21, P2211))+Math.log(Math.pow(2*P21*P22, P2212))+Math.log(Math.pow(P22*P22, P2222));
			if (LK > maxLK) {
				maxLK = LK;
				FP11 = EP11;
				FP12 = EP12;
				FP21 = EP21;
				FP22 = EP22;
			}

			//System.out.println("log likeihood " +maxLK);
		}
			//System.out.println(FP11+" "+FP12+" "+FP21+" "+FP22);
			//System.out.println(p1+" "+p2+" "+q1+" "+q2);
	}

	public void set(int genNum, int genType) {
		if (genType == 1111)
			P1111 = genNum;
		if (genType == 1112)
			P1112 = genNum;
		if (genType == 1122)
			P1122 = genNum;
		if (genType == 1211)
			P1211 = genNum;
		if (genType == 1212)
			P1212 = genNum;
		if (genType == 1222)
			P1222 = genNum;
		if (genType == 2211)
			P2211 = genNum;
		if (genType == 2212)
			P2212 = genNum;
		if (genType == 2222)
			P2222 = genNum;
	}
	public static void main (String[] args) {
		EMHaplo eh=new EMHaplo(1, 2);
		eh.set(10, 1111);
		eh.set(15, 1112);
		eh.set(5, 1122);
		eh.set(10, 1211);
		eh.set(50, 1212);
		eh.set(13, 1222);
		eh.set(3, 2211);
		eh.set(13, 2212);
		eh.set(10, 2222);
		eh.estimate();
	}

}