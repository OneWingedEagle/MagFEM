package fem;
import static java.lang.Math.PI;
import static java.lang.Math.*;
import math.Vect;
import math.util;


public class Region {

	public int dim,BHnumber,lamBNumber,currentIndex,unknownCurrentIndex,curMap1;
	private String material,name;
	private double mu0=PI*4e-7;
	private Vect mur,sigma,nu,J,M,er;
	public double phase0,freq,omega,ro,thermalCoef,deltaT,beta=0,wireRes,nloop;
	public Vect pois,yng,shear;
	private int firstElement,lastElement;
	public  boolean hasJ,stranded,circuit, hasM,isConductor,isNonLinear,deformable,MS,rotor,thermal,isotElast;
	public double Jz,windingSurf,NtS,current,currentp,currCoef1,terminalVoltage0,terminalVoltage,terminalVoltagep,inducedVoltage;
	
	public Region(int dim)
	{
		this.dim=dim;
		 mur=new Vect().ones(dim);
		 er=new Vect().ones(dim);
		 nu=mur.times(mu0).inv();
		 sigma=new Vect(3);
		 J=new Vect(dim);
		 M=new Vect(dim);
		 
	}
	
	public void setFirstEl(int first){
		this.firstElement=first;
		
	}
	
	public void setLastEl(int last){
		this.lastElement=last;
		
	}
	
	public int getFirstEl(){
		return this.firstElement;
		
	}
	
	public int getLastEl(){
		return this.lastElement;
	}
	
	public void setMaterial(String mat){
		this.material=mat;
		
	}
	public void setName(String name){
		this.name=name;
		
			
	}
	
	public String getMaterial(){
		return this.material;
		
	}
	
	public String getName(){
		return this.name;
		
	}
	public void setJ(Vect J){

		this.J=J.deepCopy();
		if(J.norm()==0)
		this.hasJ=false;
		else
		this.hasJ=true;
	}
	
	public Vect getJ(){

		return this.J.deepCopy();
	
	}

	public void setM(Vect M){
	
		this.M=M.deepCopy().times(1);
		if(M.norm()==0)
		this.hasM=false;
		else
		this.hasM=true;
	
	}
	
	public Vect getM(){

		return this.M.deepCopy();
	
	}
		
	
	public void setSigma(Vect sigma){
		this.sigma=sigma.deepCopy();
		if(sigma.norm()>0)
		this.isConductor=true;
		
	}
	
	public Vect getSigma(){

		return this.sigma.deepCopy();
	}
	
	
	public Vect getEr(){

		return this.er.deepCopy();
	}
	

	public void setMur(Vect mur){

	this.mur=mur.deepCopy();
	
//	this.mur=new Vect().ones(dim);
		this.nu=this.mur.times(this.mu0).inv();
		
	
				
	}
	
	public void setEr(Vect er){

		this.er=er.deepCopy();
							
		}
	
	public void setNu(Vect nu){

		 this.nu=nu.deepCopy();
	
	}
	
	public void setNu(double nu){

		 this.nu=new Vect(dim);
		 for(int k=0;k<dim;k++)
				this.nu.el[k]=nu;	
	
	}
	public Vect getNu(){

		return this.nu.deepCopy();
	
	}
	
	public Vect getMu(){

		return this.nu.inv();
	
	}
	
	public Vect getMur(){

		return this.mur.deepCopy();
	
	}
	
	public void setNonLinear(boolean b){
		this.isNonLinear=b;
				
	}
	
	public void setPois(Vect pois){
		this.pois=pois.deepCopy();
		
	}
	public void setYng(Vect yng){
		this.yng=yng.deepCopy();
	}
	public void setShear(Vect shear){
		this.shear=shear.deepCopy();
	}

	public void setRo(double ro){
		this.ro=ro;
	}
	
	public void setThermalCoef(double alpha){
		this.thermalCoef=alpha;
		if(alpha!=0) thermal=true;
	}
	
	public double getThermalCoef(){
		return this.thermalCoef;
	}
	
	
	public void setDeltaT(double dT){
		this.deltaT=dT;
	
	}
	
	public double getDeltaT(){
		return this.deltaT;
	
	}
	
	public Vect getPois(){
		if(this.pois==null) return new Vect(dim);;
	return	this.pois.deepCopy();
		
	}

	public Vect getYng(){
		if(this.yng==null)  return new Vect(dim);;;
		return	this.yng.deepCopy();
			
		}
	public Vect getShear(){
		if(this.shear==null) return new Vect(dim);;
		return	this.shear.deepCopy();
			
		}
	
	public void setBeta(double beta){

		 this.beta=beta;
	
	}

	public double getRo(){
		return	this.ro;
			
		}
	
/*	public void setFactJx(double fj){
		
		this.factJ=fj;

	}
	public void setFactJcosx(double dt){
	
		//this.factJ=1*signum(cos(phase0+beta+elAng)-1.1*util.saw(70*(elAng)/PI/2));
		
		//this.factJ=signum(phase0+1);
		//this.factJ=1.05*util.saw(70*elAng/PI/2);
	this.factJ=1*cos(phase0+beta+omega*dt);//-.2*sin(2*(phase0+beta+elAng))+.1*cos(4*(phase0+beta+elAng))-.05*sin(6*(phase0+beta+elAng))+.03*cos(8*(phase0+beta+elAng));
	
	}
	
	public double getFactJx(){
		return	this.factJ;

		}*/
	
/*	public void setJAt(double d){
		
		//this.factJ=1*signum(cos(phase0+beta+elAng)-1.1*util.saw(70*(elAng)/PI/2));
		
		//this.factJ=signum(phase0+1);
		//this.factJ=1.05*util.saw(70*elAng/PI/2);
	this.factJ=J0.times*cos(phase0+beta+omega*dt);//-.2*sin(2*(phase0+beta+elAng))+.1*cos(4*(phase0+beta+elAng))-.05*sin(6*(phase0+beta+elAng))+.03*cos(8*(phase0+beta+elAng));
	
	}*/
	public void setDeformable(boolean b){
		this.deformable=b;
	}
	public int getNumbElements(){
		return this.lastElement-this.firstElement+1;
	}
	
	public double getWireRes(){
		return	this.wireRes;
			
		}
	public void setWireRes(double R){
			this.wireRes=R;
			
		}
	
	public void setFreq(double f){
			this.freq=f;
			this.omega=2*PI*f;
			
		}
	
	public Region deepCopy(){
		Region rc=new Region(this.dim);
		rc.dim=this.dim;
		rc.name=this.name;
		rc.material=this.material;
		rc.firstElement=this.firstElement;
		rc.lastElement=this.lastElement;
		
		return rc;
	}
}
	
	
