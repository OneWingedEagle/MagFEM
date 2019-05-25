package meshFactory;

import math.util;

public class Geometry  {



	public double[][] blockBoundary;
	public double[][] minMeshLeft,minMeshRight,baseLeft,baseRight;
	public boolean[][] discLeft,discRight;
	public double minMesh=10,scaleFactor=1000,base=1.5;
	public int nBlocks,nBoundary;
	public String[] blockName;
	public int coordinate=0;
	public Geometry()  {}
	
	public Geometry(double[][] bb)  {
		
	
		this.blockBoundary=bb;
		this.nBoundary=bb[0].length;
		nBlocks=bb.length;
		discLeft=new boolean[nBlocks][this.nBoundary];
		discRight=new boolean[nBlocks][this.nBoundary];
		
		minMeshLeft=new double[nBlocks][this.nBoundary];
		minMeshRight=new double[nBlocks][this.nBoundary];
		
		baseLeft=new double[nBlocks][this.nBoundary];
		baseRight=new double[nBlocks][this.nBoundary];
		blockName=new String[nBlocks];
		for(int i=0;i<nBlocks;i++)
			for(int j=0;j<nBoundary;j++){
				discLeft[i][j]=true;
				discRight[i][j]=true;
				
				minMeshLeft[i][j]=minMesh;
				minMeshRight[i][j]=minMesh;
				
				baseLeft[i][j]=base;
				baseRight[i][j]=base;
				
			}
		for(int i=0;i<nBlocks;i++){
		
			blockName[i]="block"+(i+1);
		}
	

	
	}

	
}
	
	
		
	


