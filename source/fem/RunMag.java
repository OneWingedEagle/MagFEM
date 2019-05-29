package fem;
import java.io.File;

import java.text.DecimalFormat;

import femSolver.StaticElectricSolver;
import main.Main;
import math.Vect;
import math.util;

public class RunMag {

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
			int ix=0;
			
			int nBegin=model.nBegin;
			int nEnd=model.nEnd;
			

			int inc=model.nInc;


			
			String currentFolder=model.resultFolder;
			String fluxFolder=model.resultFolder;
			String dispFolder=model.resultFolder;
			
			if(model.saveForce){
		
			
					dispFolder = folder+"\\forces";
				File dfolder = new File(dispFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}

			if(model.saveFlux){
					fluxFolder = folder+"\\fluxes";
				
				File dfolder = new File(fluxFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();

			}
			util.pr(model.saveJe);
			if(model.saveJe){
				currentFolder = folder+"\\currents";
			
			File dfolder = new File(currentFolder);
			if(dfolder.exists())
				util.deleteDir(dfolder);
			dfolder.mkdir();

		}

			
			if(nTsteps>1)
				model.writeFiles=false;
			else
				model.writeFiles=true;

			boolean writeFiles=model.writeFiles;
			
			
			Vect x=new Vect();

			main.gui.lbX[0].setText("rot ang. ");
			main.gui.lbX[1].setText("Trq. ");

			
			model.solveCoils();
			

				for(int step=nBegin;step<=nEnd;step+=inc){
			
			
				
				main.gui.tfX[0].setText((step)+"/"+nEnd);
				
				model.currentTimeStep=step;
				
				model.setJ0(model.getCurrentTime());	
				
				if(model.phiCoils!=null)
				model.phiCoils[0].current=Math.cos(2*model.freq*step*model.dt);

					if(!model.loadFlux){
						
						
							model.setMagBC();
							
						
							if(step==nBegin){
								model.writer.reportData(model);
							}

							
							if(!model.nonLin /*|| i==nBegin*//* ||model.Q!=null*/){
								
							
									if(!model.loadPrevMag){
										model.saveAp();		
										x=model.solveMagLin(step-nBegin,x);	
									

									}
									else
									{
							
									}

								}
								else
									x=this.xp;
						
								
							if(model.nonLin ){

								if(x==null) x=new Vect(model.numberOfUnknowns);
									if(model.analysisMode>0 && step!=nBegin )	{	
								
										model.saveAp();				
									}
						
									x=model.solveNonLinear(x,true,step-nBegin);
								}

						if(model.analysisMode>0)
							model.setJe();
				
						if(x!=null)
							this.xp=x.deepCopy();

	/*						if(model.saveFlux){
								String fluxFile = fluxFolder+"\\fluxes\\flux"+step+".txt";
								if(step==nBegin)
									model.writeMesh(fluxFolder+"\\fluxes\\bun"+step+".txt");
			
								model.writeB(fluxFile);
								
							//	model.writer.writeA_as_flux(model,fluxFolder+"\\fluxA"+i+".txt");
								
		
							}*/


					//	model.resetReluctForce();
					//	model.setReluctForce();
						//
						//model.setMSForce();
						
						
				}


					
			
					/*if(model.dim==2)
						model.setTorque(0,model.rm,1);
						else
							model.setTorque(model.r1,1,1);

					
					
					util.pr("torque >>>>>>>"+model.TrqZ);*/
					
				//	T.el[ix++]=model.TrqZ;
					
								//************************************************
				//	writeFiles=true;

					if(model.phiCoils==null)
						model.writeJ0(folder+"\\J"+step+".txt");
					

					model.setSolution(x);	
					
					double loss=0;
					model.setJe();
					
					loss=model.writer.outputLoss(model,model.resultFolder+"\\outputs.txt",step,model.getCurrentTime());

					if(model.saveJe){
						String JeFile =  currentFolder+"\\Je"+step+".txt";
	
	
						model.writeJe(JeFile);
					
						if(step==nBegin)
							model.writeMesh(currentFolder+"\\bun"+step+".txt");
						
		
					}
					
					if(model.saveFlux){

						String fluxFile = fluxFolder+"\\flux"+step+".txt";
						if(step==nBegin)
							model.writeMesh(fluxFolder+"\\bun"+step+".txt");
						model.writeB(fluxFile);
						
					}
					
				
					
					model.writer.outputEnergies(model,model.resultFolder+"\\outputs.txt",step,model.getCurrentTime());

					
					T.el[ix++]=loss;
				
					
					if(model.saveForce){
					//	String dispFile =dispFolder+"\\force"+mangs+".txt";
						String forceFile =model.resultFolder+"\\force"+step+".txt";
						model.writeNodalField(forceFile,model.forceCalcMode);
					}


					main.gui.tfX[1].setText(this.formatter.format(model.TrqZ));

					if(model.solver.terminate){

						break;
					}

			if(model.solver.terminate) break;


			}
			//	util.pr(model.main.gui.iccgArea.getText());

				//util.plot(T);
			//	T.show();
				
				Vect errs=new Vect(model.solver.totalIter);
			
				for(int i=0;i<errs.length;i++)
					errs.el[i]=model.solver.errs.get(i);

				//errs.show();

				util.plot("error",errs.el,"ICCG Convergence");

			//	util.plot(errs);




	}

}
