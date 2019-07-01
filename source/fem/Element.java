package fem;
import math.Mat;
import math.Vect;
import math.util;


public class Element {
	private int nRegion;
	public byte dim;
	private int[] vertexNumb;
	private int[] edgeNumb;
	private boolean[] edgeReverse;
	public int[] edgeXYNumb;
	private Vect nu,sigma,B,J,Je,M;
	private boolean hasJ,hasM,isConductor,deformable,nonlin,MS,thermal;
	

	public Element(String type){

		if(type.equals("triangle") ){
			this.vertexNumb=new int[3];
			this.edgeNumb=new int[3];
			edgeReverse=new boolean[3];
			dim=2;
		}
			
		else if(type.equals("quadrangle") ){
			this.vertexNumb=new int[4];
			this.edgeNumb=new int[4];
			edgeReverse=new boolean[4];
			dim=2;
		}
		else if(type.equals("tetrahedron") ){
			this.vertexNumb=new int[4];
			this.edgeNumb=new int[6];
			edgeReverse=new boolean[6];
			dim=3;
		}
		
		else if(type.equals("prism") ){
			this.vertexNumb=new int[6];
			this.edgeNumb=new int[9];
			edgeReverse=new boolean[9];
			dim=3;
		}
		else if(type.equals("hexahedron") ){
			this.vertexNumb=new int[8];
			this.edgeNumb=new int[12];
			edgeReverse=new boolean[12];

			dim=3;
		}
		
		else if(type.equals("pyramid") ){
			this.vertexNumb=new int[5];
			this.edgeNumb=new int[8];
			edgeReverse=new boolean[8];
			dim=3;
		}
		
		this.B=new Vect(dim);
		this.nu=new Vect().ones(dim);
		
	}
	
	public void setJ(Vect J){
		this.J=J.deepCopy();
		this.hasJ=true;

	}
	public Vect getJ(){
		return this.J.deepCopy();

	}

	public void setM(Vect M){
		this.M=M.deepCopy();
		this.hasM=true;

	}
	
	public void setNu(Vect nu){

			this.nu=nu.deepCopy();
	}

	public Vect getM(){
		if(hasM)
		return this.M.deepCopy();
		else return new Vect(dim);

	}


	public void setB(Vect B){

		this.B=B;

	}
	
	public void setB(int k,double Bk){

		this.B.el[k]=Bk;

	}
	

	public Vect getB(){
		return this.B.deepCopy();

	}

	public void setJe(Vect Je){
		this.Je=Je;

	}

	public Vect getJe(){
	
		if(this.Je==null) return new Vect(3);
		return this.Je.deepCopy();

	}


	public Vect getNu(){
		return this.nu.deepCopy();

	}
	

	public void setRegion(int nr){
		this.nRegion=nr;

	}
	
	public int getRegion(){
		
		return this.nRegion;
	}

	public void setSigma(Vect sigma){

		if(sigma.norm()>0) {
			this.sigma=sigma.deepCopy();
			this.isConductor=true;
			Je=new Vect(3);
		}
		
		else
			this.isConductor=false;


	}


	public Vect getSigma(){
		
		if(sigma!=null)
		return this.sigma.deepCopy();
		else
			return new Vect(3);
	}
	
	public double getSigmaZ(){
		
		if(isConductor)
		return this.sigma.el[2];
		else  return 0.0;

	}

	
	
	public int getDim(){
		
		return this.dim;
	}
	
	public void setEdgeNumb(int[] ne){
		int nEdge=ne.length;
		edgeNumb=new int[nEdge];
		for(int i=0;i<nEdge;i++){
			edgeNumb[i]=ne[i];
		}
	}
	
	public int[] getEdgeNumb(){
		int nEdge=edgeNumb.length;
		int[] ne=new int[nEdge];
		for(int i=0;i<nEdge;i++)
			ne[i]=edgeNumb[i];
		return  ne;
	}
	
	
	public void setEdgeReverse(boolean[] er){
		int nEdge=er.length;
		edgeReverse=new boolean[nEdge];
		for(int i=0;i<nEdge;i++){
			edgeReverse[i]=er[i];
		}
	}
	
	public boolean[] getEdgeReverse(){
		int nEdge=edgeReverse.length;
		boolean[] er=new boolean[nEdge];
		for(int i=0;i<nEdge;i++)
			er[i]=edgeReverse[i];
		return  er;
	}
	public void setVertNumb(int[] nv){
		int nVert=nv.length;
		 vertexNumb=new int[nVert];
		for(int i=0;i<nVert;i++)
			vertexNumb[i]=nv[i];
	}
	
	public int[] getVertNumb(){
		int nVert=vertexNumb.length;
		int[] nv=new int[nVert];
		for(int i=0;i<nVert;i++)
			nv[i]=vertexNumb[i];
		return  nv;
	}
	
	public void setEdgeNumb(int j, int ne){

			edgeNumb[j]=ne;
	}
	
	public int getEdgeNumb(int j){

		return  edgeNumb[j];
	}
	
	public void setVertNumb(int j, int nv){

		vertexNumb[j]=nv;
	}
	
	public int getVertNumb(int j){

		return  vertexNumb[j];
	}
	
	public void setHasJ(boolean b){
		this.hasJ=b;
	}
	public void setHasM(boolean b){
		this.hasM=b;
	}
	
	public void setHasThermal(boolean b){
		this.thermal=b;
	}
	
	public boolean  isThermal(){
		return this.thermal;
	}

	public void setNonlin(boolean b){
		this.nonlin=b;
	}
	
	public boolean hasJ(){
		return this.hasJ;
	}
	
	public boolean hasM(){
		return this.hasM;
	}
	public boolean hasMS(){
		return this.MS;
	}
	public boolean isDeformable(){
		return this.deformable;
	}
	
	public boolean isConductor(){
		return (this.isConductor);
	}
	
	public void setConductor(boolean b){
		this.isConductor=b;
	}
	
	public boolean isNonlin(){
		return this.nonlin;
	}
	

	
}
