package math;

import static java.lang.Math.*;
import io.Loader;

import java.awt.Color;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.Frame;
import java.io.File;
import java.text.DecimalFormat;
import java.util.Random;

import javax.swing.JFileChooser;
import javax.swing.JFrame;










import materialData.CurrentWaveForm;

import org.math.plot.Plot2DPanel;

		
		public class CLN {
			static DecimalFormat df=new DecimalFormat("0.00000E00");
			static String regex="[ : ,=\\t]+";

			public CLN(){}
			
			public static void main(String[] args) throws Exception{
				
				CLN  cln=  new CLN();
	
				cln.ladder_R_netowrk();

			}
	

	
	private static double ladder_R_netowrk() {
	

	int nStages = 5;

	Vect Rs=new Vect(nStages);
	Vect Rp=new Vect(nStages);
	
	Rs.el[0]=2;
	Rp.el[0]=1;
	
	for (int k = 1; k <nStages; k++){
		
		if((k-1)%3==0){
			int d=(k-1)/3;
			Rs.el[k]=(4*d+2);
		}else
			Rs.el[k]=1.;
		
		if((k-2)%3==0){
			int d=(k-2)/3;
			Rp.el[k]=1./(4*d+4);
		}else
			Rp.el[k]=1.;
	}
	
	Rs.hshow();
	Rp.hshow();

	double Req;
	int last = nStages - 1;
	Req = Rs.el[last];
	Req+= Rp.el[last];

	for (int k = nStages - 2; k >= 0; k--) {
		double G=0;
		if (Rp.el[k]>0)
			G += 1. / Rp.el[k];
		if (Req>0)
			G += 1. / Req;
	
		Req = Rs.el[k];
		if (G>0)
			Req += 1./G;

	}
	
	double myE=Req;
	
	util.pr(myE);
	util.pr(Math.E);
	
	util.pr("erro(%)= "+100*abs(Math.E-myE)/Math.E);
	return Req;

	}
}
