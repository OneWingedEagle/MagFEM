package fem;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.awt.Color;
import java.io.File;

import java.text.DecimalFormat;

import femSolver.StaticElectricSolver;
import main.Main;
import math.Vect;
import math.VectComp;
import math.util;

public class RunMagAC {

	private DecimalFormat formatter=new DecimalFormat("0.00");

	private Vect xp;
	int[] nn=new int[1];

	public static void main(String[] args){
		new Main();
	}

	public void runMag(Model model, Main main){
		

			String folder=model.resultFolder;


			int nTsteps=model.nTsteps;
		
			Vect T=new Vect(nTsteps);
			Vect freqs=new Vect(nTsteps);
			int ix=0;
			
			int nBegin=model.nBegin;
			int nEnd=model.nEnd;
			

			int inc=model.nInc;

			if(model.AC) nEnd=nBegin;

			
			String dispFolder=model.resultFolder;
			String fluxFolder=model.resultFolder;

			if(model.saveForce){
			
/*			
					dispFolder = folder+"\\forces";
				File dfolder = new File(dispFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();*/

			}

			if(model.saveFlux){
	/*				fluxFolder = folder+"\\fluxes";
				
				File dfolder = new File(fluxFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();*/

			}


			
			if(nTsteps>1)
				model.writeFiles=false;
			else
				model.writeFiles=true;

			boolean writeFiles=model.writeFiles;
			
			
			VectComp xc=new VectComp();

			main.gui.lbX[0].setText("rot ang. ");
			main.gui.lbX[1].setText("Trq. ");


		
			model.solveCoils();



			for(int i=nBegin;i<=nEnd;i+=inc){
					
				double t0=0;
		
		
				main.gui.tfX[2].setText((i)+"/"+nEnd);

				model.setJ0(t0);	
				
				if(ix>0){
				model.freq/=Math.pow(10, .5);
				freqs.el[ix]=model.freq;
				
				}
					
					model.setMagBC();
							
							if(i==nBegin){
								model.writer.reportData(model);
							}

							xc=model.femSolver.solveMagAC(model, i-nBegin);	
				

							int m=xc.length;
							Vect vr1=new Vect(m);
							Vect vm1=new Vect(m);

							Vect vr2=new Vect(m);
							Vect vm2=new Vect(m);
							for(int k=0;k<m;k++){
								if(k<model.numberOfUnknownEdges){
									vr1.el[k]=xc.el[k].re;
									vm1.el[k]=xc.el[k].im;
									
									vr2.el[k]=-xc.el[k].im*2*PI;
									vm2.el[k]=xc.el[k].re*2*PI;
								}
								else{
									vr2.el[k]=xc.el[k].re;
									vm2.el[k]=xc.el[k].im;
								}
								
							}
							
							
							double lossRe=0;
							double lossIm=0;

								if(model.saveFlux){
									model.setSolution(vr1);	

									String fluxFile = folder+"\\fluxRe"+i+".txt";
									model.writeB(fluxFile);
								}
								

								model.setSolution(vr2);	
								model.setJe();
								
								if(model.saveJe){
									String JeFile =  folder+"\\JeRe"+i+".txt";
	
				
									model.writeJe(JeFile);
						
								}
										
	
								for(int ir=1;ir<=model.numberOfRegions;ir++){
									if(model.region[ir].isConductor){
								lossRe=model.obtainLoss(ir);
								util.pr("Real Joule Loss of region "+ir+" [Re]=  "+lossRe);
									}
								}

								if(model.saveFlux){
									model.setSolution(vm1);	
									String fluxFile = folder+"\\fluxIm"+i+".txt";
									model.writeB(fluxFile);
								}
								
								model.setSolution(vm2);	
								model.setJe();	
								if(model.saveJe){
										String JeFile =  folder+"\\JeIm"+i+".txt";
						
										model.writeJe(JeFile);
						
								}
								
								for(int ir=1;ir<=model.numberOfRegions;ir++){
									if(model.region[ir].isConductor){
										lossIm=model.obtainLoss(ir);
								util.pr("Joule Loss of region "+ir+" [Imag]=  "+lossIm);
									}
								}
							
			
						T.el[ix++]=(lossRe+lossIm)/2;



					main.gui.tfX[1].setText(this.formatter.format(model.TrqZ));

					if(model.solver.terminate) break;

				}
			
			//util.plot(freqs,T);
			//T.show();
		
				Vect errs=new Vect(model.solver.totalIter);

				for(int i=0;i<errs.length;i++)
					errs.el[i]=model.solver.errs.get(i);
				
				util.plot("error",errs.el,"ICCG Convergence");
			//	util.plot(errs);

	

	}

}
